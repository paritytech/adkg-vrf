#![cfg_attr(not(feature = "std"), no_std)]

use crate::bls::threshold::ThresholdVk;
use crate::dkg::aggregator::TranscriptAggregator;
use crate::pvss::SecretSharing;
use crate::sig_agg::SignatureAggregator;
use ark_ec::pairing::Pairing;
use dkg::transcript;

pub mod bls;
pub mod dkg;
pub mod koe;
/// Threshold Verifiable Unpredictable Function (VUF) scheme.
/// Produces an unpredictable output by aggregating a threshold number of vanilla BLS signatures on the input.
///
/// The scheme comprises 2 parts:
/// 1. a Distributed Key Generation (DKG) protocol that produces some data for a set of BLS signers, and
/// 2. a BLS signature aggregation scheme that leverages the data produced by the DKG
///    to aggregate the signatures from a subset of the signers into a threshold signature,
///    and additionally produce a VUF output.
///
/// An interesting property of the scheme is that the signers are not required to participate
/// in the protocol in any way other than producing vanilla BLS signatures.
/// That allows to transform any deployed BLS signature scheme, where the same message is being signed by multiple signers,
/// into a threshold scheme or a randomness beacon.
///
/// The implementation follows the notes by Alistair Stewart:
/// 1. https://hackmd.io/3968Gr5hSSmef-nptg2GRw
/// 2. https://hackmd.io/xqYBrigYQwyKM_0Sn5Xf4w
/// TODO: is there a paper?

/// Aggregatable Publicly Verifiable Secret Sharing Scheme
// mod old_dkg;
pub mod pvss;
pub mod sig_agg;
pub mod straus;
pub mod utils;
mod bs_dkg;

mod hash_to_curve;
pub use hash_to_curve::PairingWithG1Map;

/// Verified aggregated secret shared to a list of signers with a specified threshold.
/// Has all the data required to aggregate or verify threshold signatures for a single committee.
/// Contains signers' BLS public keys only in G2,
/// so isn't immediately suitable for the evolving committee case.
#[derive(Clone, Debug)]
pub struct VerifiedSharing<C: Pairing> {
    /// Key material.
    secret_sharing: SecretSharing<C>,
    /// The list of signers and the threshold.
    params: pvss::Params<C>,
}

impl<C: Pairing> VerifiedSharing<C> {
    pub fn config(&self) -> pvss::Config<C> {
        self.params.config.clone()
    }

    pub fn verification_key(&self) -> ThresholdVk<C> {
        ThresholdVk::from_share(&self.secret_sharing)
    }

    pub fn aggregation_key(&self) -> SignatureAggregator<C> {
        SignatureAggregator::new(
            &self.params.signer_pks,
            self.secret_sharing.bgpk.clone(),
            self.params.config.clone(),
        )
    }

    pub fn into_keys(self) -> (ThresholdVk<C>, SignatureAggregator<C>) {
        let vk = self.verification_key();
        let agg = self.aggregation_key();
        (vk, agg)
    }
}


pub type BlsDkg = dkg::Dkg<ark_bls12_381::Bls12_381>;
pub type BlsSignerPk = ark_bls12_381::G2Affine;
pub type BlsTranscriptAggregator = TranscriptAggregator<ark_bls12_381::Bls12_381>;
pub type DkgTranscript = transcript::Transcript<ark_bls12_381::Bls12_381>;

// must have
// TODO: Fiat-Shamir
// TODO: cofactors/subgroup checks
// TODO: return results
// TODO: ark-substrate

// nice to have
// TODO: CP proofs
// TODO: integration test
// TODO: half-aggregation
// TODO: bench for logn = 16, 20

// nice to consider
// TODO: IBE
// TODO: resharing?
// TODO: backsharing
// TODO: multiple Cs?

// TODO: test single signer, t = 1
// TODO: test t = n
// TODO: test multiple dealings
#[cfg(test)]
mod tests {
    use crate::bls::threshold::ThresholdVk;
    use crate::bls::vanilla::BlsSigner;
    use crate::dkg::transcript::Transcript;
    use crate::dkg::Dkg;
    use crate::pvss;

    use crate::sig_agg::SignatureAggregator;
    use ark_bls12_381::{Bls12_381, G1Affine};
    use ark_ec::pairing::Pairing;
    use ark_ec::{AffineRepr, CurveGroup};
    use ark_std::rand::Rng;
    use ark_std::test_rng;
    use ark_std::vec::Vec;

    // Returns threshold verification and aggregation keys
    pub fn simulate_pvss<C: Pairing, R: Rng>(
        signers_pks: Vec<C::G2Affine>,
        t: usize,
        rng: &mut R,
    ) -> (ThresholdVk<C>, SignatureAggregator<C>) {
        let pvss = pvss::Params::<C>::new(signers_pks.clone(), t).unwrap();
        let share = pvss.deal(rng).unwrap();
        let tvk = ThresholdVk::from_share(&share.payload);
        let bgpk = share.payload.bgpk;
        let sig_aggregator =
            SignatureAggregator::<C>::new(signers_pks.as_slice(), bgpk, pvss.config);
        (tvk, sig_aggregator)
    }

    #[test]
    fn it_works() {
        let rng = &mut test_rng();

        let (n, t) = (7, 5);
        let signers: Vec<BlsSigner<Bls12_381>> = (0..n).map(|_| BlsSigner::new(rng)).collect();
        let signers_pks: Vec<_> = signers.iter().map(|s| s.bls_pk_g2).collect();

        let dealers: Vec<_> = (0..3).map(|_| BlsSigner::<Bls12_381>::new(rng)).collect();
        let dealer_pks: Vec<G1Affine> = dealers.iter().map(|d| d.bls_pk_g1).collect();

        let dkg =
            Dkg::<Bls12_381>::new(signers_pks.clone(), t, dealer_pks.clone(), dealer_pks.len())
                .unwrap();

        let transcripts: Vec<Transcript<Bls12_381>> = dealers
            .into_iter()
            .map(|dealer| {
                dkg.deal_and_sign(rng, (dealer.sk, dealer.bls_pk_g1))
                    .unwrap()
            })
            .collect();

        assert!(dkg.verify(&transcripts[0], rng).is_ok());

        let agg_transcript = Dkg::<Bls12_381>::aggregate(transcripts);

        assert!(dkg.verify(&agg_transcript, rng).is_ok());

        let keys = dkg.finalize(agg_transcript, rng).unwrap();

        let config = keys.config();
        let threshold_vk = ThresholdVk::from_share(&keys.secret_sharing);
        let sig_aggregator =
            SignatureAggregator::<Bls12_381>::new(&signers_pks, keys.secret_sharing.bgpk, config);

        let message = BlsSigner::<Bls12_381>::hash_to_g1("message".as_bytes()).into_group();
        let sigs: Vec<_> = signers.iter().map(|s| s.sign_g1(message)).collect();

        let threshold_sig_n =
            sig_aggregator.check_then_aggregate(message.into_affine(), sigs.clone());
        let vuf_n = threshold_vk.vuf_unoptimized(&threshold_sig_n, message);

        let sigs_t: Vec<_> = sigs.into_iter().take(t).collect();
        let threshold_sig_t =
            sig_aggregator.check_then_aggregate(message.into_affine(), sigs_t.clone());
        let vuf_t = threshold_vk.vuf_unoptimized(&threshold_sig_t, message);

        assert_eq!(vuf_n, vuf_t);

        let (ss, epk) = threshold_vk.initiate_key_exchange(b"message", rng);
        let ss_ = epk.complete_key_exchange(&threshold_sig_t);
        assert_eq!(ss, ss_);
        let ss_ = epk.complete_key_exchange(&threshold_sig_t);
        assert_eq!(ss, ss_);
    }
}
