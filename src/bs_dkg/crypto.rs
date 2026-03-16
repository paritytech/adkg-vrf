use crate::bls::vanilla::StandaloneSig;
use crate::pvss;
use crate::utils::BarycentricDomain;
use ark_ec::pairing::Pairing;
use ark_ec::{AffineRepr, CurveGroup, VariableBaseMSM};
use ark_ff::Zero;

pub struct EvolvingCommitteeTpk<C: Pairing> {
    pub c: C::G1Affine,
    pub g1: C::G1Affine,
    pub g2: C::G2Affine,
}

pub struct EcAggThresholdSig<C: Pairing> {
    pub asig: C::G1Affine,
    pub apk_g1: C::G1Affine,
    pub apk_g2: C::G2Affine,
    pub abgpk_tweaked: C::G2Affine,
}

impl<C: Pairing> EvolvingCommitteeTpk<C> {
    pub fn with_c(c: C::G1Affine) -> Self {
        Self {
            c,
            g1: C::G1Affine::generator(),
            g2: C::G2Affine::generator(),
        }
    }

    pub fn verify_sig(&self, sig: &EcAggThresholdSig<C>, msg: C::G1Affine, h2_pred: C::G2Affine) {
        // BLS aggregate public keys consistency across `G1` and `G2`.
        // `apk_g1 = ask.g1` and `apk_g2 = ask.g2` for some `ask`.
        // `e(g1, apk_g2) = e(apk_g1, g2)`
        assert_eq!(
            C::pairing(self.g1, sig.apk_g2),
            C::pairing(sig.apk_g1, self.g2)
        );
        // BLS signature verification against a public key in `G2`.
        // `e(asig, g2) = e(msg, apk_g2)`
        assert_eq!(
            C::pairing(sig.asig, self.g2),
            C::pairing(msg, sig.apk_g2)
        );
        // Aggregation consistency.
        // Let `h2 = sh.g2`, then `sh.pk_j = sh.(sk_j.g2) = sk_j.(sh.g2) = sk_j.h2`.
        // `bgpk_j = f(w^j).g2 + sh.pk_j = f(w^j).g2 + sk_j.h2`
        // `tweak_j = sk_j.(h2_pred - h2) = sk_j.h2_pred - sk_j.h2`
        // `bgpk_tweaked_j = bgpk_j + tweak_j = f(w^j).g2 + sk_j.h2_pred`
        // `e(g1, abgpk_tweaked) = e(C, g2) + e(apk_g1, h2_pred)`
        assert_eq!(
            C::pairing(self.g1, sig.abgpk_tweaked),
            C::multi_pairing(&[self.c, sig.apk_g1], &[self.g2, h2_pred])
        );
        // TODO: e(C, g2) is constant
        // TODO: prepare h2_pred
    }
}

pub fn aggregate_sigs<C: Pairing>(
    bls_sigs: Vec<StandaloneSig<C>>,
    pks_g1: Vec<C::G1Affine>,
    bgpks: Vec<C::G2Affine>,
    config: &pvss::Config<C>,
) -> EcAggThresholdSig<C> {
    let lis = BarycentricDomain::of_size(config.domain, config.n)
        .lagrange_basis_at(C::ScalarField::zero());
    let sigs: Vec<_> = bls_sigs.iter().map(|s| s.sig).collect();
    let pks_g2: Vec<_> = bls_sigs.iter().map(|s| s.pk).collect();
    let asig = C::G1::msm(&sigs, &lis).unwrap().into_affine();
    let apk_g1 = C::G1::msm(&pks_g1, &lis).unwrap().into_affine();
    let apk_g2 = C::G2::msm(&pks_g2, &lis).unwrap().into_affine();
    let abgpk = C::G2::msm(&bgpks, &lis).unwrap().into_affine();
    EcAggThresholdSig {
        asig,
        apk_g1,
        apk_g2,
        abgpk_tweaked: abgpk,
    }
}