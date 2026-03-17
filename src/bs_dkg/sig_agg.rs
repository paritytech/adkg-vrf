use crate::bls::vanilla::StandaloneSig;
use crate::bs_dkg::crypto::EcAggThresholdSig;
use crate::pvss;
use crate::utils::BarycentricDomain;
use ark_ec::pairing::Pairing;
use ark_ec::{CurveGroup, VariableBaseMSM};
use ark_ff::Zero;
use ark_poly::EvaluationDomain;
use hashbrown::HashMap;

#[derive(Clone)]
struct AugmentedEcSig<C: Pairing> {
    pub sig: C::G1Affine,
    pub pk_g1: C::G1Affine,
    pub pk_g2: C::G2Affine,
    pub bgpk_tweaked: C::G2Affine,
}

pub fn aggregate_ec_sigs<C: Pairing>(
    augmented_sigs: Vec<Option<AugmentedEcSig<C>>>,
    config: &pvss::Config<C>,
) -> EcAggThresholdSig<C> {
    assert_eq!(augmented_sigs.len(), config.n);
    let mut bitmask: Vec<bool> = augmented_sigs.iter().map(|o| o.is_some()).collect();
    bitmask.resize(config.domain.size(), false);
    let set_bits_count = bitmask.iter().filter(|b| **b).count();
    assert!(set_bits_count >= config.t);
    let lis = BarycentricDomain::from_subset(config.domain, &bitmask)
        .lagrange_basis_at(C::ScalarField::zero());
    let augmented_sigs: Vec<AugmentedEcSig<C>> = augmented_sigs.into_iter().flatten().collect();
    let sigs: Vec<C::G1Affine> = augmented_sigs.iter().map(|s| s.sig).collect();
    let pks_g1: Vec<C::G1Affine> = augmented_sigs.iter().map(|s| s.pk_g1).collect();
    let pks_g2: Vec<C::G2Affine> = augmented_sigs.iter().map(|s| s.pk_g2).collect();
    let bgpks: Vec<C::G2Affine> = augmented_sigs.iter().map(|s| s.bgpk_tweaked).collect();
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
struct PkExt<C: Pairing> {
    j: usize,
    pk_g1: C::G1Affine,
    bgpk_tweaked: C::G2Affine,
}

pub struct EcSigAgg<C: Pairing> {
    /// PVSS config
    config: pvss::Config<C>,
    /// map pk_g2_j -> (j, pk_g1_j, bgpk_tweaked_j)
    pk_ext: HashMap<C::G2Affine, PkExt<C>>,
}

impl<C: Pairing> EcSigAgg<C> {
    pub fn new(
        pks_g2: Vec<C::G2Affine>,
        pks_g1: Vec<C::G1Affine>,
        bgpks_tweaked: Vec<Option<C::G2Affine>>,
        config: pvss::Config<C>,
    ) -> Self {
        let pk_ext: HashMap<_, _> = pks_g2.into_iter()
            .enumerate()
            .zip(pks_g1.into_iter())
            .zip(bgpks_tweaked.into_iter())
            .filter_map(|(((j, pk_g2), pk_g1), bgpk)|
                bgpk.map(|bgpk_tweaked| (pk_g2, PkExt { j, pk_g1, bgpk_tweaked })))
            .collect();
        Self {
            config,
            pk_ext,
        }
    }

    pub fn aggregate(&self, sigs: Vec<StandaloneSig<C>>) -> EcAggThresholdSig<C> {
        let mut augmented_sigs = vec![None; self.config.n];
        sigs.into_iter().for_each(|sig| {
            let StandaloneSig { sig, pk } = sig;
            if let Some((j, sig_ext)) = self.pk_ext.get(&pk).map(|pk_ext| {
                let sig_ext = AugmentedEcSig {
                    sig,
                    pk_g1: pk_ext.pk_g1,
                    pk_g2: pk,
                    bgpk_tweaked: pk_ext.bgpk_tweaked,
                };
                (pk_ext.j, sig_ext)
            }) {
                augmented_sigs[j] = Some(sig_ext);
            }
        });
        aggregate_ec_sigs(augmented_sigs, &self.config)
    }
}