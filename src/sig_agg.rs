use ark_ec::pairing::Pairing;
use ark_ec::VariableBaseMSM;
use ark_poly::EvaluationDomain;
use ark_std::Zero;
use ark_std::{vec, vec::Vec};
use hashbrown::HashMap;
use std::error::Error;
use std::fmt::Display;
use std::iter;

use crate::bls::threshold::AggThresholdSig;
use crate::bls::vanilla::BlsSig;
use crate::pvss;
use crate::utils::BarycentricDomain;
use ark_ec::CurveGroup;

#[derive(Debug)]
pub struct InterpolationError {
    n_evals: usize,
    degree: usize,
}
impl Error for InterpolationError {}

impl Display for InterpolationError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        self.fmt(f)
    }
}

/// Evaluations at `0` of Lagrange basis polynomials of the set of "indices" of the points.
///
/// Let `n = config.n` and `S` be the set of roots of unity corresponding to the indices of the `points`.
/// `S = (w^{i_k})` where `1 <= i_1,...,i_s <= n` are such that `points[i_k].is_some(), k=1,...,s`,
/// and `(L_1,...,L_s) = L_S` -- the Lagrange basis polynomials of `S`: `L_j(w^{i_k}) = 1 iff j = k, j = 1,...,s`.
///
/// Returns `(L_1(0),...,L_s(0)` and the corresponding points `(points[i_1], ..., points[i_s])`.
/// Requires `config.t <= s <= config.n`.
///
/// If `points[i] = p(w^{i})` for a degree `t < s` polynomial `p`, then `p(0) = L_1(0).points[i_1] + ... + L_s(0).points[i_s]`.
/// As an optimization, `s = config.t`.
pub fn evaluate_lagrange_basis_at_0<T, C: Pairing>(opt_points: Vec<Option<T>>, config: &pvss::Config<C>) -> Result<(Vec<C::ScalarField>, Vec<T>), InterpolationError> {
    debug_assert_eq!(opt_points.len(), config.n);
    let n_evals = opt_points.iter().flatten().count();
    if n_evals < config.t {
        return Err(InterpolationError{ n_evals, degree: config.t - 1 });
    }
    let mut set_bits = 0;
    let bitmask: Vec<bool> = opt_points.iter()
        .scan(&mut set_bits, |counter, opt| {
            if **counter == config.t {
                return None;
            }
            if opt.is_some() {
                **counter += 1;
            }
            Some(opt.is_some())
        })
        .chain(iter::repeat(false))
        .take(config.domain.size())
        .collect();
    debug_assert!(set_bits == config.t);
    let lis = BarycentricDomain::from_subset(config.domain, &bitmask)
        .lagrange_basis_at(C::ScalarField::zero());
    let points: Vec<T> = opt_points.into_iter()
        .flatten()
        .take(config.t)
        .collect();
    debug_assert_eq!(lis.len(), config.t);
    debug_assert_eq!(points.len(), config.t);
    Ok((lis, points))
}

/// Assuming `evals[i] = p(w^{i}).g2`, returns `p(0).g2`
pub fn evaluate_at_0_in_g2<C: Pairing>(evals_opt: Vec<Option<C::G2>>, config: &pvss::Config<C>) -> Result<C::G2, InterpolationError> {
    let (lis_at_zero, evals) = evaluate_lagrange_basis_at_0(evals_opt, config)?;
    let evals = C::G2::normalize_batch(&evals);
    let res = C::G2::msm(&evals, &lis_at_zero).unwrap();
    Ok(res)
}

/// To aggregate vanilla BLS signatures, they have to be:
/// 1. equipped with the signers' `bgpk`s
/// 2. arranged in the correct PVSS order.
/// That means that the aggregator has to know the signers' BLS public keys in order,
/// and the corresponding `bgpk`s

/// Aggregates a (not less than) threshold amount of signatures with `bgpk`s.
pub fn aggregate_augmented_sigs<C: Pairing>(
    augmented_sigs: Vec<Option<AggThresholdSig<C>>>,
    config: &pvss::Config<C>,
) -> Option<AggThresholdSig<C>> {
    let (lis, augmented_sigs) = evaluate_lagrange_basis_at_0(augmented_sigs, config).ok()?;
    let bls_sigs: Vec<_> = augmented_sigs
        .iter()
        .map(|s| s.bls_sig_with_pk.sig)
        .collect();
    let bls_pks: Vec<_> = augmented_sigs
        .iter()
        .map(|s| s.bls_sig_with_pk.pk)
        .collect();
    let bgpks: Vec<_> = augmented_sigs.iter().map(|s| s.bgpk).collect();
    let asig = C::G1::msm(&bls_sigs, &lis).unwrap().into_affine();
    let apk = C::G2::msm(&bls_pks, &lis).unwrap().into_affine();
    let abgpk = C::G2::msm(&bgpks, &lis).unwrap().into_affine();
    Some(AggThresholdSig {
        bls_sig_with_pk: BlsSig { sig: asig, pk: apk },
        bgpk: abgpk,
    })
}


/// Converts vanilla BLS signatures to the (threshold) aggregatable counterparts.
/// In order to do that knows the mapping of the vanilla BLS public keys
/// to the corresponding `bgpk` and the index in the signers list.
pub struct SignatureAggregator<C: Pairing> {
    // PVSS config
    pub(crate) config: pvss::Config<C>,
    // map bls_pk_j -> (bgpk_j, j)
    pub(crate) pks_mapping: HashMap<C::G2Affine, (C::G2Affine, usize)>,
}

impl<C: Pairing> SignatureAggregator<C> {
    /// BLS public keys and the `bgpk`s in the PVSS order.
    pub fn new(
        signer_pks: &[C::G2Affine],
        bgpks: Vec<C::G2Affine>,
        config: pvss::Config<C>,
    ) -> Self {
        let pks_mapping: HashMap<_, _> = signer_pks
            .iter()
            .cloned()
            .zip(bgpks)
            .enumerate()
            .map(|(j, (signer_pk_j, bgpk_j))| (signer_pk_j, (bgpk_j, j)))
            .collect();
        Self {
            config,
            pks_mapping,
        }
    }

    /// If the public key is recognized, returns the signer's index `j` and adds `bgpk_j` to the signature.
    pub fn augment_sig(&self, sig: BlsSig<C>) -> Option<(usize, AggThresholdSig<C>)> {
        self.pks_mapping.get(&sig.pk).map(|(bgpk, j)| (*j, AggThresholdSig {
            bls_sig_with_pk: sig,
            bgpk: *bgpk,
        }))
    }

    /// Aggregates signatures. Checks that there is a threshold of signers, but doesn't verify the individual signatures.
    /// If a public key is not recognized, the signature is dropped. Signatures from
    pub fn aggregate_wo_checking(&self, sigs: Vec<BlsSig<C>>) -> AggThresholdSig<C> {
        let mut augmented_sigs = vec![None; self.pks_mapping.len()];
        sigs.into_iter().for_each(|sig| {
            let (j, s) = self.augment_sig(sig).unwrap();
            augmented_sigs[j] = Some(s);
        });
        aggregate_augmented_sigs(augmented_sigs, &self.config).unwrap()
    }

    /// Checks that each signature verifies and comes from a legit signer.
    /// Checks that the threshold is met. TODO: Doesn't allow duplicate signatures?
    pub fn check_then_aggregate(&self, message: C::G1Affine, sigs: Vec<BlsSig<C>>) -> AggThresholdSig<C> {
        let mut session = self.start_session(message);
        session.append_verify_sigs(sigs.clone());
        let augmented_sigs = session.finalize();
        let threshold_sig = aggregate_augmented_sigs(augmented_sigs, &self.config).unwrap();
        threshold_sig
    }

    pub fn start_session(&self, message: C::G1Affine) -> Session<C> {
        Session {
            g2: self.config.g2.into_affine(),
            message,
            pks: &self.pks_mapping,
            augmented_sigs: vec![None; self.pks_mapping.len()],
        }
    }

    pub fn aggregate_tweaks(&self, _message: C::G2Affine, sigs: &[(C::G2Affine, C::G2Affine)]) -> Vec<Option<C::G2Affine>> {
        let mut tweaked_bgpks = vec![None; self.config.n];
        for (sig, pk) in sigs {
            // TODO: check the sig
            let (bgpk, j) = *self.pks_mapping.get(pk).unwrap();
            tweaked_bgpks[j] = Some((bgpk + sig).into_affine()); // TODO: batch
        }
        tweaked_bgpks
    }
}

pub struct Session<'a, C: Pairing> {
    // to verify BLS sigs with the keys in G2
    g2: C::G2Affine,
    // the message on that signatures are being aggregated
    message: C::G1Affine,
    // map bls_pk_j -> (bgpk_j, j)
    pks: &'a HashMap<C::G2Affine, (C::G2Affine, usize)>,
    /// `Vec` of length `n` accumulating the augmented signatures stored at the right index.
    augmented_sigs: Vec<Option<AggThresholdSig<C>>>,
}

impl<'a, C: Pairing> Session<'a, C> {
    pub fn finalize(self) -> Vec<Option<AggThresholdSig<C>>> {
        self.augmented_sigs
        // params.aggregate_augmented_sigs(self.augmented_sigs)
    }

    /// Signatures MUST
    /// 1. valid on the message
    /// 2. from a known pk
    /// duplicates allowed
    /// TODO: return result of indices
    pub fn append_verify_sig(&mut self, sig: BlsSig<C>) {
        let (bgpk, j) = {
            let x = self.pks.get(&sig.pk);
            assert!(x.is_some());
            x.unwrap().clone()
        };
        assert!(self.augmented_sigs[j].is_none());
        sig.verify_unoptimized(self.message.into(), self.g2);
        self.augmented_sigs[j] = Some(AggThresholdSig {
            bls_sig_with_pk: sig,
            bgpk,
        })
    }

    pub fn append_verify_sigs(&mut self, sigs: Vec<BlsSig<C>>) {
        sigs.into_iter().for_each(|s| self.append_verify_sig(s));
    }
}
