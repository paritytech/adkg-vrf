use crate::bls::vanilla::BlsSig;
use crate::bs_dkg::crypto::{aggregate_ec_sigs, EvolvingCommitteeAggSig, EvolvingCommitteeSig};
use crate::pvss;
use ark_ec::pairing::Pairing;
use ark_ec::CurveGroup;
use hashbrown::HashMap;


struct PkExt<C: Pairing> {
    j: usize,
    pk_g1: C::G1Affine,
    bgpk: C::G2Affine,
}

pub struct EcSigAgg<C: Pairing> {
    /// Epoch id
    sid: u64,
    /// PVSS config
    config: pvss::Config<C>,
    /// map pk_g2_j -> (j, pk_g1_j, bgpk_tweaked_j)
    pk_ext: HashMap<C::G2Affine, PkExt<C>>,
    bgpk_delta: C::G2,
}

impl<C: Pairing> EcSigAgg<C> {
    pub fn new(
        sid: u64,
        pks_g2: Vec<C::G2Affine>,
        pks_g1: Vec<C::G1Affine>,
        bgpks: Vec<Option<C::G2Affine>>,
        bgpk_delta: C::G2,
        config: pvss::Config<C>,
    ) -> Self {
        let pk_ext: HashMap<_, _> = pks_g2.into_iter()
            .enumerate()
            .zip(pks_g1.into_iter())
            .zip(bgpks.into_iter())
            .filter_map(|(((j, pk_g2), pk_g1), bgpk)|
                bgpk.map(|bgpk_tweaked| (pk_g2, PkExt { j, pk_g1, bgpk: bgpk_tweaked })))
            .collect();
        Self {
            sid,
            config,
            pk_ext,
            bgpk_delta,
        }
    }

    pub fn aggregate(&self, sigs: Vec<BlsSig<C>>) -> Option<EvolvingCommitteeAggSig<C>> {
        (sigs.len() >= self.config.t).then(|| ())?;
        let mut augmented_sigs = vec![None; self.config.n];
        sigs.into_iter().for_each(|sig| {
            if let Some((j, sig_ext)) = self.pk_ext.get(&sig.pk).map(|pk_ext| {
                let sig_ext = EvolvingCommitteeSig {
                    sig: sig.sig,
                    pk_g1: pk_ext.pk_g1,
                    pk_g2: sig.pk,
                    bgpk: pk_ext.bgpk,
                };
                (pk_ext.j, sig_ext)
            }) {
                augmented_sigs[j] = Some(sig_ext);
            }
        });
        let mut asig = aggregate_ec_sigs(augmented_sigs, &self.config)?;
        asig.bgpk = (asig.bgpk - self.bgpk_delta).into_affine();
        Some(EvolvingCommitteeAggSig {
            sid: self.sid,
            asig,
        })
    }
}