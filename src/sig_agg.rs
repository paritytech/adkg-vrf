use ark_ec::pairing::Pairing;
use ark_ec::VariableBaseMSM;
use ark_poly::EvaluationDomain;
use ark_std::Zero;
use ark_std::{vec, vec::Vec};
use hashbrown::HashMap;

use crate::bls::threshold::AggThresholdSig;
use crate::bls::vanilla::StandaloneSig;
use crate::pvss;
use crate::utils::BarycentricDomain;
use ark_ec::CurveGroup;


/// To aggregate vanilla BLS signatures, they have to be:
/// 1. equipped with the signers' `bgpk`s
/// 2. arranged in the correct PVSS order.
/// That means that the aggregator has to know the signers' BLS public keys in order,
/// and the corresponding `bgpk`s

/// Aggregates a (not less than) threshold amount of signatures with `bgpk`s.
pub fn aggregate_augmented_sigs<C: Pairing>(
    augmented_sigs: Vec<Option<AggThresholdSig<C>>>,
    config: &pvss::Config<C>,
) -> AggThresholdSig<C> {
    assert_eq!(augmented_sigs.len(), config.n);
    let mut bitmask: Vec<bool> = augmented_sigs.iter().map(|o| o.is_some()).collect();
    bitmask.resize(config.domain.size(), false);
    let set_bits_count = bitmask.iter().filter(|b| **b).count();
    assert!(set_bits_count >= config.t);
    let lis = BarycentricDomain::from_subset(config.domain, &bitmask)
        .lagrange_basis_at(C::ScalarField::zero());
    let augmented_sigs: Vec<AggThresholdSig<C>> =
        augmented_sigs.into_iter().flatten().collect();
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
    AggThresholdSig {
        bls_sig_with_pk: StandaloneSig { sig: asig, pk: apk },
        bgpk: abgpk,
    }
}


/// Converts vanilla BLS signatures to the threshold aggregatable counterparts.
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
    pub fn start_session(&self, message: C::G1Affine) -> Session<C> {
        Session {
            g2: self.config.g2.into_affine(),
            message,
            pks: &self.pks_mapping,
            augmented_sigs: vec![None; self.pks_mapping.len()],
        }
    }

    pub fn aggregate(&self, message: C::G1Affine, sigs: Vec<StandaloneSig<C>>) -> AggThresholdSig<C> {
        let mut session = self.start_session(message);
        session.append_verify_sigs(sigs.clone());
        let augmented_sigs = session.finalize();
        let threshold_sig = aggregate_augmented_sigs(augmented_sigs, &self.config);
        threshold_sig
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
    pub fn append_verify_sig(&mut self, sig: StandaloneSig<C>) {
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

    pub fn append_verify_sigs(&mut self, sigs: Vec<StandaloneSig<C>>) {
        sigs.into_iter().for_each(|s| self.append_verify_sig(s));
    }
}
