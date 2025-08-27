use ark_ec::pairing::Pairing;
use ark_ec::{CurveGroup, PrimeGroup, VariableBaseMSM};
use ark_ff::Zero;
use ark_poly::EvaluationDomain;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::vec::Vec;
use hashbrown::HashMap;

use crate::bls::threshold::AggThresholdSig;
use crate::bls::vanilla::StandaloneSig;
use crate::old_dkg::verifier::TranscriptVerifier;
use crate::utils::BarycentricDomain;

/// Aggregatable Publicly Verifiable Secret Sharing scheme (aPVSS) for sharing a secret key `f(0).g1` in G1,
/// corresponding to the public key `f(0).g2` in G2.
///
/// There are 2 types of participants:
/// 1. a fixed list of signers, identified with their BLS public keys in G2, and
/// 2. any number of dealers, whose authentication is the problem of a higher-level protocol.
///
/// A dealer samples a secret and produces a transcript containing shares of the secret, each encrypted to the corresponding signer,
/// together with a publicly verifiable proof of validity of the ciphertexts.
/// Transcripts with contributions from different dealers can be aggregated into a single verifiable transcript.
/// The scheme is secure (vaguely, that means the parameters produced are secure),
/// if the final aggregated transcript is valid, and contains a contribution from a single honest dealer.
///
/// *A fun property* of the scheme is that signers don't have to use (or even decrypt) their shares.
/// Instead, anyone can blindly use the ciphertexts to produce proofs that the threshold number of signers have signed.

pub mod dealer;
pub mod verifier;
pub mod transcript;

//TODO: move bls_pks out?
/// Parameters of an aPVSS instantiation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Ceremony<'a, C: Pairing, D: EvaluationDomain<C::ScalarField>> {
    /// The number of signers.
    pub n: usize,
    /// The threshold, i.e. the minimal number of signers required to reconstruct the shared secret.
    pub t: usize,
    /// The signers' bls public keys in G2.
    /// **It's critical that proofs of possession are checked for these keys.**
    pub bls_pks: &'a [C::G2Affine],
    /// An FFT-friendly multiplicative subgroup of the field of size not less than `n`.
    /// The evaluation points are the first `n` elements of the subgroup: `x_j = w^j, j = 0,...,n-1`,
    /// where `w` is the generator of the subgroup.
    pub domain: D,
    /// Generator of G1.
    pub g1: C::G1,
    /// Generator of G2.
    pub g2: C::G2,
}

/// Useful data produced by the protocol:
/// - encrypted shares of the secret key,
/// - the corresponding threshold public key, and
/// - a pair of points with the same discrete logarithm.
///
/// The secret key being shared among the signers is `f(0).g2` for some polynomial `f`.
/// `f(0).g1` is the public key, corresponding to the shared secret key. The share of the signer `j` is `f(w^j).g2`.
/// `(h1, h2)` are points in G1xG2 with the same discrete logarithm, i.e. `h1 = sh.g1` and `h2 = sh.g2` for some `sh`.
/// `bgpk_j = f(w^j).g2 + sh.pk_j, j = 0,...,n-1`.
/// Then `(bgpk_j, h2)` is the ElGamal encryption of the point `f(w^j).g2` with `pk_j` for the ephemeral secret `sh`.
#[derive(Clone, Debug, PartialEq, Eq, CanonicalDeserialize, CanonicalSerialize)]
//TODO: check visibility
//TODO: better name
pub struct DkgResult<C: Pairing> {
    /// The public key corresponding to the shared secret key.
    /// `c = f(0).g1`
    pub(crate) c: C::G1Affine,
    /// Shares of the secret, encrypted to the signers.
    /// `bgpk_j = f(w^j).g2 + sh.pk_j, j = 0,...,n-1`
    pub bgpk: Vec<C::G2Affine>,
    /// `h1 = sh.g1`
    pub(crate) h1: C::G1Affine,
    /// `h2 = sh.g2`
    pub(crate) h2: C::G2Affine,
}

impl<'a, C: Pairing, D: EvaluationDomain<C::ScalarField>> Ceremony<'a, C, D> {
    pub fn setup(t: usize, bls_pks: &'a [C::G2Affine]) -> Self {
        let n = bls_pks.len();
        assert!(t <= n);
        //todo: test t = 1, t = n
        Self {
            n,
            t,
            bls_pks,
            domain: D::new(n).unwrap(),
            g1: C::G1::generator(),
            g2: C::G2::generator(),
        }
    }

    pub fn verifier(&self) -> TranscriptVerifier<C> {
        TranscriptVerifier::new_with_domain(self.domain, self.n, self.t)
    }

    // TODO: args are not any more aggregatable
    pub fn aggregate_augmented_sigs(&self, augmented_sigs: Vec<Option<AggThresholdSig<C>>>) -> AggThresholdSig<C> {
        assert_eq!(augmented_sigs.len(), self.n);
        let mut bitmask: Vec<bool> = augmented_sigs.iter().map(|o| o.is_some()).collect();
        bitmask.resize(self.domain.size(), false);
        let set_bits_count = bitmask.iter().filter(|b| **b).count();
        assert!(set_bits_count >= self.t);
        let lis = BarycentricDomain::from_subset(self.domain, &bitmask)
            .lagrange_basis_at(C::ScalarField::zero());
        let augmented_sigs: Vec<AggThresholdSig<C>> = augmented_sigs.into_iter()
            .flatten()
            .collect();
        let bls_sigs: Vec<_> = augmented_sigs.iter().map(|s| s.bls_sig_with_pk.sig).collect();
        let bls_pks: Vec<_> = augmented_sigs.iter().map(|s| s.bls_sig_with_pk.pk).collect();
        let bgpks: Vec<_> = augmented_sigs.iter().map(|s| s.bgpk).collect();
        let asig = C::G1::msm(&bls_sigs, &lis).unwrap().into_affine();
        let apk = C::G2::msm(&bls_pks, &lis).unwrap().into_affine();
        let abgpk = C::G2::msm(&bgpks, &lis).unwrap().into_affine();
        AggThresholdSig {
            bls_sig_with_pk: StandaloneSig { sig: asig, pk: apk },
            bgpk: abgpk,
        }
    }

    pub fn aggregator(&self, final_share: DkgResult<C>) -> crate::agg::SignatureAggregator<C> {
        let pks: HashMap<_, _> = self.bls_pks.iter()
            .cloned()
            .zip(final_share.bgpk)
            .enumerate()
            .map(|(j, (bls_pk_j, bgpk_j))| (bls_pk_j, (bgpk_j, j)))
            .collect();
        crate::agg::SignatureAggregator {
            g2: self.g2.into_affine(),
            pks,
        }
    }
}

impl<C: Pairing> DkgResult<C> {
    pub fn merge_with(self, mut others: Vec<Self>) -> Self {
        others.push(self);
        Self::merge(&others)
    }

    pub fn merge(keys: &[Self]) -> Self {
        let n = keys[0].bgpk.len();
        Self {
            c: (keys.iter().map(|key| key.c).sum::<C::G1>()).into_affine(),
            // TODO: affine conversions
            bgpk: (0..n).map(|j| {
                keys.iter()
                    .map(|key| key.bgpk[j])
                    .sum::<C::G2>()
                    .into_affine()
            }).collect(),
            h1: keys.iter().map(|key| key.h1).sum::<C::G1>().into_affine(),
            h2: keys.iter().map(|key| key.h2).sum::<C::G2>().into_affine(),
        }
    }
}
