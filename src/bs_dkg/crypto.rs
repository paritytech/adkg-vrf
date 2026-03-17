use crate::bs_dkg::BsDkg;
use crate::hash_to_curve::CurveWithPairingAndHash;
use crate::pvss;
use crate::utils::BarycentricDomain;
use ark_ec::pairing::Pairing;
use ark_ec::{AffineRepr, CurveGroup, VariableBaseMSM};
use ark_ff::Zero;
use ark_poly::EvaluationDomain;

pub struct EvolvingCommitteePk<C: Pairing> {
    pub c: C::G1Affine,
    pub g1: C::G1Affine,
    pub g2: C::G2Affine,
}

#[derive(Clone, Debug)]
pub(crate) struct EvolvingCommitteeSig<C: Pairing> {
    pub(crate) sig: C::G1Affine,
    pub(crate) pk_g1: C::G1Affine,
    pub(crate) pk_g2: C::G2Affine,
    pub(crate) bgpk: C::G2Affine,
}

pub struct EvolvingCommitteeAggSig<C: Pairing> {
    pub(crate) sid: u64,
    pub(crate) asig: EvolvingCommitteeSig<C>,
}

impl<C: CurveWithPairingAndHash> EvolvingCommitteePk<C> {
    pub fn with_c(c: C::G1Affine) -> Self {
        Self {
            c,
            g1: C::G1Affine::generator(),
            g2: C::G2Affine::generator(),
        }
    }

    pub fn verify(&self, sig: &EvolvingCommitteeAggSig<C>, msg: &[u8]) {
        let msg_in_g1 = C::hash_to_g1(msg).unwrap();
        let h2_pred = BsDkg::<C>::h2_of(sig.sid);
        self.verify_point(sig, msg_in_g1, h2_pred)
    }

    pub fn verify_point(&self, sig: &EvolvingCommitteeAggSig<C>, msg: C::G1Affine, h2_pred: C::G2Affine) {
        let sig = &sig.asig;
        // BLS aggregate public keys consistency across `G1` and `G2`.
        // `apk_g1 = ask.g1` and `apk_g2 = ask.g2` for some `ask`.
        // `e(g1, apk_g2) = e(apk_g1, g2)`
        assert_eq!(
            C::pairing(self.g1, sig.pk_g2),
            C::pairing(sig.pk_g1, self.g2)
        );
        // BLS signature verification against a public key in `G2`.
        // `e(asig, g2) = e(msg, apk_g2)`
        assert_eq!(
            C::pairing(sig.sig, self.g2),
            C::pairing(msg, sig.pk_g2)
        );
        // Aggregation consistency.
        // Let `h2 = sh.g2`, then `sh.pk_j = sh.(sk_j.g2) = sk_j.(sh.g2) = sk_j.h2`.
        // `bgpk_j = f(w^j).g2 + sh.pk_j = f(w^j).g2 + sk_j.h2`
        // `tweak_j = sk_j.(h2_pred - h2) = sk_j.h2_pred - sk_j.h2`
        // `bgpk_tweaked_j = bgpk_j + tweak_j = f(w^j).g2 + sk_j.h2_pred`
        // `e(g1, abgpk_tweaked) = e(C, g2) + e(apk_g1, h2_pred)`
        assert_eq!(
            C::pairing(self.g1, sig.bgpk),
            C::multi_pairing(&[self.c, sig.pk_g1], &[self.g2, h2_pred])
        );
        // TODO: e(C, g2) is constant
        // TODO: prepare h2_pred
    }
}

pub fn aggregate_ec_sigs<C: Pairing>(
    augmented_sigs: Vec<Option<EvolvingCommitteeSig<C>>>,
    config: &pvss::Config<C>,
) -> EvolvingCommitteeSig<C> {
    assert_eq!(augmented_sigs.len(), config.n);
    let mut bitmask: Vec<bool> = augmented_sigs.iter().map(|o| o.is_some()).collect();
    bitmask.resize(config.domain.size(), false);
    let set_bits_count = bitmask.iter().filter(|b| **b).count();
    assert!(set_bits_count >= config.t);
    let lis = BarycentricDomain::from_subset(config.domain, &bitmask)
        .lagrange_basis_at(C::ScalarField::zero());
    let augmented_sigs: Vec<EvolvingCommitteeSig<C>> = augmented_sigs.into_iter().flatten().collect();
    let sigs: Vec<C::G1Affine> = augmented_sigs.iter().map(|s| s.sig).collect();
    let pks_g1: Vec<C::G1Affine> = augmented_sigs.iter().map(|s| s.pk_g1).collect();
    let pks_g2: Vec<C::G2Affine> = augmented_sigs.iter().map(|s| s.pk_g2).collect();
    let bgpks: Vec<C::G2Affine> = augmented_sigs.iter().map(|s| s.bgpk).collect();
    let asig = C::G1::msm(&sigs, &lis).unwrap().into_affine();
    let apk_g1 = C::G1::msm(&pks_g1, &lis).unwrap().into_affine();
    let apk_g2 = C::G2::msm(&pks_g2, &lis).unwrap().into_affine();
    let abgpk = C::G2::msm(&bgpks, &lis).unwrap().into_affine();
    EvolvingCommitteeSig {
        sig: asig,
        pk_g1: apk_g1,
        pk_g2: apk_g2,
        bgpk: abgpk,
    }
}