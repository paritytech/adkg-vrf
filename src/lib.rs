#![cfg_attr(not(feature = "std"), no_std)]

use crate::bls::threshold::ThresholdVk;
use crate::dkg::aggregator::TranscriptAggregator;
use crate::pvss::SecretSharing;
use crate::sig_agg::SignatureAggregator;
use crate::utils::BarycentricDomain;
use ark_ec::pairing::Pairing;
use ark_ec::CurveGroup;
use ark_ec::VariableBaseMSM;
use ark_std::Zero;
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

#[derive(Clone, Debug)]
pub struct VerifiedSharingWithG1Keys<C: Pairing> {
    verified_sharing: VerifiedSharing<C>,
    signers_g1: Vec<C::G1Affine>,
    h2_pred: C::G2Affine,
}

impl<C: Pairing> VerifiedSharingWithG1Keys<C> {

    pub fn tweak_msg(&self) -> C::G2 {
        self.h2_pred - self.verified_sharing.secret_sharing.h2
    }

    pub fn tweak_bgpks(&self, sigs: &[(C::G2Affine, C::G2Affine)]) -> Vec<Option<C::G2Affine>> {
        let tweak_msg = self.tweak_msg().into_affine();
        let tweaked_bgpks = self.verified_sharing.aggregation_key().aggregate_tweaks(tweak_msg, sigs);
        tweaked_bgpks
    }

    pub fn tweak(self, sigs: &[(C::G2Affine, C::G2Affine)]) -> TweakedSharing<C> {
        // TODO: verify
        let tweaked_bgpks = self.tweak_bgpks(sigs);
        TweakedSharing {
            secret_sharing: self.verified_sharing,
            signers_g1: self.signers_g1,
            tweaked_bgpks,
            h2_pred: self.h2_pred,
        }
    }
}

/// Verified aggregated (related) secrets shared to `2` consequent committees of signers.
#[derive(Clone, Debug)]
pub struct VerifiedSharingAndBack<C: Pairing> {
    back_sharing: VerifiedSharingWithG1Keys<C>,
    next_sharing: VerifiedSharingWithG1Keys<C>,
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

pub struct TweakedSharing<C: Pairing> {
    /// Secret shared to signers with keys in G2.
    secret_sharing: VerifiedSharing<C>,
    /// Signers' keys in G1 in the same order as `self.secret_sharing.params.signer_pks`
    signers_g1: Vec<C::G1Affine>,
    tweaked_bgpks: Vec<Option<C::G2Affine>>,
    h2_pred: C::G2Affine,
}

impl<C: Pairing> TweakedSharing<C> {
    /// Computes the threshold public key delta. This is (f0(0) - f1_back(0)).g1
    /// `deg(f0) = deg(f1_back) = t0 - 1`.
    /// Thus both the current sharing and the back-shared one should have `t0` tweaked bgpks at the same positions.
    /// TODO: check for that
    pub fn compute_delta(&self, bs_next: &Self) -> C::G2Affine {
        let bgpk_deltas: Vec<C::G2> = self.tweaked_bgpks.iter()
            .zip(bs_next.tweaked_bgpks.iter())
            .map(|(curr, bs_next)| bs_next.unwrap() - curr.unwrap())
            .collect();
        let bgpk_deltas = C::G2::normalize_batch(&bgpk_deltas);
        let bs_next_config = &bs_next.secret_sharing.params.config;
        let lis_at_zero = BarycentricDomain::of_size(bs_next_config.domain, bs_next_config.n)
            .lagrange_basis_at(C::ScalarField::zero());
        C::G2::msm(&bgpk_deltas, &lis_at_zero).unwrap()
            .into_affine()
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

    use crate::bs_dkg::crypto::{aggregate_sigs, EcAggThresholdSig, EvolvingCommitteeTpk};
    use crate::bs_dkg::{BsDkg, Committee};
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

    // TODO:
    // 1. signers' pks in g1
    // 2. signer's produce tweaks
    // 3. tweaks are verified, aggregated and applied
    // 4. new key material type
    // 5. new signature type
    // 6. aggregation for that
    // basically 2 cryptosuites, ideally with a fall-back option

    #[test]
    fn test_back_sharing() {
        let rng = &mut test_rng();

        // Protocol and the participants
        let (n, t) = (7, 5);

        let signers_0: Vec<BlsSigner<Bls12_381>> = (0..n).map(|_| BlsSigner::new(rng)).collect();
        let signers_pks_g2_0: Vec<_> = signers_0.iter().map(|s| s.bls_pk_g2).collect();
        let signers_pks_g1_0: Vec<_> = signers_0.iter().map(|s| s.bls_pk_g1).collect();

        let signers_1: Vec<BlsSigner<Bls12_381>> = (0..n).map(|_| BlsSigner::new(rng)).collect();
        let signers_pks_g2_1: Vec<_> = signers_1.iter().map(|s| s.bls_pk_g2).collect();
        let signers_pks_g1_1: Vec<_> = signers_1.iter().map(|s| s.bls_pk_g1).collect();

        let first = Committee {
            params: pvss::Params::<Bls12_381>::new(signers_pks_g2_0.clone(), t).unwrap(),
            signers_g1: signers_pks_g1_0.clone(),
        };
        let next = Committee {
            params: pvss::Params::<Bls12_381>::new(signers_pks_g2_1.clone(), t).unwrap(),
            signers_g1: signers_pks_g1_1.clone(),
        };
        let bs_dkg = BsDkg::start(first);

        // Deals secret shares to the epoch #1 committee (no-one to backshare to)
        let transcript = bs_dkg.deal_first(signers_0[0].clone(), rng).unwrap();
        let ss_0 = bs_dkg.verify_first(transcript, rng);
        let config_0 = ss_0.verified_sharing.params.config.clone();
        let h2_pred_0 = ss_0.h2_pred;
        let (fc_tpk_0, sig_aggregator_0) = ss_0.clone().verified_sharing.into_keys();
        let ec_tpk = EvolvingCommitteeTpk::with_c(fc_tpk_0.c);

        // Tests a threshold signature at epoch 0
        let msg = BlsSigner::<Bls12_381>::hash_to_g1("message".as_bytes()).into_group();
        let sigs: Vec<_> = signers_0.iter().map(|s| s.sign_g1(msg)).collect();
        let agg_sig_0 = sig_aggregator_0.aggregate_wo_checking(sigs.clone());
        fc_tpk_0.verify_unoptimized(&agg_sig_0, msg);

        // Tweaks the key material (`bgpks` and `h2`)
        let tweak_msg_0 = ss_0.tweak_msg();
        let tweaks_0: Vec<_> = signers_0.iter()
            .map(|s| (s.sign_g2(tweak_msg_0), s.bls_pk_g2))
            .collect();
        let ss_tweaked_0 = ss_0.tweak(&tweaks_0);

        // TODO: write the aggregator for the evolving committee scheme
        let bgpk_mod_0: Vec<_> = ss_tweaked_0.tweaked_bgpks.iter().map(|x| x.unwrap()).collect();
        let agg_sig_mod_0 = aggregate_sigs(
            sigs,
            signers_pks_g1_0,
            bgpk_mod_0.clone(),
            &config_0,
        );

        ec_tpk.verify_sig(&agg_sig_mod_0, msg.into_affine(), h2_pred_0);


        let bs_dkg = bs_dkg.next(next);
        let bs_transcript = bs_dkg.deal(signers_0[0].clone(), rng).unwrap();
        let verified_bs = bs_dkg.verify(bs_transcript, rng);
        let ss_1 = verified_bs.next_sharing;
        let mut ss_1_back = verified_bs.back_sharing;
        let fc_tpk_1 = ThresholdVk::from_share(&ss_1.verified_sharing.secret_sharing);
        let h2_pred_1 = ss_1.h2_pred;

        // TWEAKS
        let tweak_msg_1 = ss_1.tweak_msg();
        let tweaks_1: Vec<_> = signers_1.iter()
            .map(|s| (s.sign_g2(tweak_msg_1), s.bls_pk_g2))
            .collect();
        let tweak_msg_back_1 = ss_1_back.tweak_msg();
        let tweaks_back_1: Vec<_> = signers_0.iter()
            .map(|s| (s.sign_g2(tweak_msg_back_1), s.bls_pk_g2))
            .collect();

        let ss_mod_1 = ss_1.tweak(&tweaks_1);
        let ss_back_mod_1 = ss_1_back.tweak(&tweaks_back_1);
        let bgpk_delta = ss_tweaked_0.compute_delta(&ss_back_mod_1);

        let bgpk_mod_1: Vec<_> = ss_mod_1.tweaked_bgpks.iter().map(|x| x.unwrap()).collect();
        let ec_tpk_2 = EvolvingCommitteeTpk::with_c(fc_tpk_1.c);
        let sigs: Vec<_> = signers_1.iter().map(|s| s.sign_g1(msg)).collect();
        let agg_sig_mod_1 = aggregate_sigs(sigs, signers_pks_g1_1, bgpk_mod_1, &config_0);

        ec_tpk_2.verify_sig(&agg_sig_mod_1, msg.into_affine(), h2_pred_1);

        let sig = EcAggThresholdSig {
            asig: agg_sig_mod_1.asig,
            apk_g1: agg_sig_mod_1.apk_g1,
            apk_g2: agg_sig_mod_1.apk_g2,
            abgpk_tweaked: (agg_sig_mod_1.abgpk_tweaked - bgpk_delta).into_affine(),
        };
        ec_tpk.verify_sig(&sig, msg.into_affine(), h2_pred_1);
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
