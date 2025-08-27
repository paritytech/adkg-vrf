use ark_ec::pairing::Pairing;
use ark_std::{vec, vec::Vec};
use hashbrown::HashMap;

use crate::bls::threshold::AggThresholdSig;
use crate::bls::vanilla::StandaloneSig;

pub struct SignatureAggregator<C: Pairing> {
    // to verify BLS sigs with the keys in G2
    pub(crate) g2: C::G2Affine,
    // map bls_pk_j -> (bgpk_j, j)
    pub(crate) pks: HashMap<C::G2Affine, (C::G2Affine, usize)>,
}

impl<C: Pairing> SignatureAggregator<C> {
    pub fn start_session(&self, message: C::G1Affine) -> Session<C> {
        Session {
            g2: self.g2,
            message,
            pks: &self.pks,
            augmented_sigs: vec![None; self.pks.len()],
        }
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
        sigs.into_iter()
            .for_each(|s| { self.append_verify_sig(s) });
    }
}

