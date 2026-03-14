#![cfg_attr(not(feature = "std"), no_std)]

use crate::bls::threshold::ThresholdVk;
use crate::bls::vanilla::StandaloneSig;
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

/// Verified aggregated secret shared to a list of signers with a specified threshold.
/// Contains all the data required to aggregate or verify threshold signatures.
pub struct VerifiedSharing<C: Pairing> {
    /// Key material.
    secret_sharing: SecretSharing<C>,
    /// The list of signers and the threshold.
    params: pvss::Params<C>,
}

/// Verified aggregated (related) secrets shared to `2` consequent committees of signers.
pub struct VerifiedBackSharing<C: Pairing> {
    back_sharing: VerifiedSharing<C>,
    next_sharing: VerifiedSharing<C>,
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

    // TODO: remove h2_pred
    pub fn tweak_h2(self, sigs: &[(C::G2Affine, C::G2Affine)], h2_pred: C::G2Affine) -> TweakedSharing<C> {
        let tweak_message = (h2_pred - self.secret_sharing.h2).into_affine();
        let tweaked_bgpks = self.aggregation_key().aggregate_tweaks(tweak_message, sigs);
        TweakedSharing {
            secret_sharing: self,
            tweaked_bgpks,
            h2_pred,
        }
    }
}

pub struct TweakedSharing<C: Pairing> {
    /// Secret shared to signers with keys in G2.
    secret_sharing: VerifiedSharing<C>,
    // /// Signers' keys in G1 in the same order as `self.secret_sharing.params.signer_pks`
    // signer_pks_g1: Vec<C::G1Affine>,
    tweaked_bgpks: Vec<Option<C::G2Affine>>,
    h2_pred: C::G2Affine,
}

impl<C: Pairing> TweakedSharing<C> {
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

pub struct Tpk<C: Pairing> {
    pub c: C::G1Affine,
    pub h2: C::G2Affine,
    pub g1: C::G1,
    pub g2: C::G2,
}

pub struct Sig<C: Pairing> {
    pub asig: C::G1Affine,
    pub apk_g1: C::G1Affine,
    pub apk_g2: C::G2Affine,
    pub abgpk: C::G2Affine,
}

impl<C: Pairing> Tpk<C> {
    pub fn verify_sig(&self, sig: &Sig<C>, m: C::G1, h2_pred: C::G2Affine) {
        assert_eq!(
            C::pairing(self.g1, sig.apk_g2),
            C::pairing(sig.apk_g1, self.g2)
        );
        assert_eq!(
            C::pairing(sig.asig, self.g2),
            C::pairing(m.into(), sig.apk_g2)
        );
        assert_eq!(
            C::pairing(self.g1.into(), sig.abgpk),
            C::multi_pairing(&[self.c, sig.apk_g1], &[self.g2.into(), h2_pred])
        );
    }
}

pub fn aggregate_sigs<C: Pairing>(
    bls_sigs: Vec<StandaloneSig<C>>,
    pks_g1: Vec<C::G1Affine>,
    bgpks: Vec<C::G2Affine>,
    config: &pvss::Config<C>,
) -> Sig<C> {
    let lis = BarycentricDomain::of_size(config.domain, config.n)
        .lagrange_basis_at(C::ScalarField::zero());
    let sigs: Vec<_> = bls_sigs.iter().map(|s| s.sig).collect();
    let pks_g2: Vec<_> = bls_sigs.iter().map(|s| s.pk).collect();
    let asig = C::G1::msm(&sigs, &lis).unwrap().into_affine();
    let apk_g1 = C::G1::msm(&pks_g1, &lis).unwrap().into_affine();
    let apk_g2 = C::G2::msm(&pks_g2, &lis).unwrap().into_affine();
    let abgpk = C::G2::msm(&bgpks, &lis).unwrap().into_affine();
    Sig {
        asig,
        apk_g1,
        apk_g2,
        abgpk,
    }
}

#[cfg(test)]
mod tests {
    use crate::bls::threshold::ThresholdVk;
    use crate::bls::vanilla::BlsSigner;
    use crate::dkg::transcript::Transcript;
    use crate::dkg::Dkg;
    use crate::{aggregate_sigs, pvss, Sig, Tpk};

    use crate::bs_dkg::BsDkg;
    use crate::sig_agg::SignatureAggregator;
    use ark_bls12_381::{Bls12_381, G1Affine, G2Affine};
    use ark_ec::pairing::Pairing;
    use ark_ec::{AffineRepr, CurveGroup};
    use ark_std::rand::Rng;
    use ark_std::test_rng;
    use ark_std::vec::Vec;
    use ark_std::UniformRand;

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

        let dealers: Vec<_> = (0..1).map(|_| BlsSigner::<Bls12_381>::new(rng)).collect();
        let dealer_pks: Vec<G1Affine> = dealers.iter().map(|d| d.bls_pk_g1).collect();

        let (n, t) = (7, 5);

        let signers_0: Vec<BlsSigner<Bls12_381>> = (0..n).map(|_| BlsSigner::new(rng)).collect();
        let signers_pks_g2_0: Vec<_> = signers_0.iter().map(|s| s.bls_pk_g2).collect();
        let signers_pks_g1_0: Vec<_> = signers_0.iter().map(|s| s.bls_pk_g1).collect();

        let dkg_0 = Dkg::<Bls12_381>::new(
            signers_pks_g2_0.clone(),
            t,
            dealer_pks.clone(),
            dealer_pks.len(),
        ).unwrap();

        let ss_0 = dkg_0.deal_and_sign(rng, (dealers[0].sk, dealers[0].bls_pk_g1)).unwrap();
        let ss_0 = dkg_0.finalize(ss_0, rng).unwrap();
        let tpk_0 = ThresholdVk::from_share(&ss_0.secret_sharing);
        let config = ss_0.params.config.clone();
        let sig_aggregator_0 = SignatureAggregator::<Bls12_381>::new(
            &signers_pks_g2_0,
            ss_0.secret_sharing.bgpk.clone(),
            config.clone(),
        );

        let message = BlsSigner::<Bls12_381>::hash_to_g1("message".as_bytes()).into_group();
        let sigs: Vec<_> = signers_0.iter().map(|s| s.sign_g1(message)).collect();
        let agg_sig_0 = sig_aggregator_0.aggregate_wo_checking(sigs.clone());
        tpk_0.verify_unoptimized(&agg_sig_0, message);

        let c_perm = tpk_0.c;
        let h2_pred_0 = G2Affine::rand(rng); // hash_to_curve(C||0)
        let tweaks_0: Vec<_> = signers_0.iter()
            .map(|s| (s.sign_g2(h2_pred_0 - ss_0.secret_sharing.h2), s.bls_pk_g2))
            .collect();
        let ss_mod_0 = ss_0.tweak_h2(&tweaks_0, h2_pred_0);
        let bgpk_mod_0: Vec<_> = ss_mod_0.tweaked_bgpks.iter().map(|x|x.unwrap()).collect();
        let tpk_mod_0 = Tpk {
            c: c_perm,
            h2: h2_pred_0,
            g1: tpk_0.g1,
            g2: tpk_0.g2,
        };
        let agg_sig_mod_0 = aggregate_sigs(
            sigs,
            signers_pks_g1_0,
            bgpk_mod_0.clone(),
            &config,
        );
        tpk_mod_0.verify_sig(&agg_sig_mod_0, message, h2_pred_0);

        let signers_1: Vec<BlsSigner<Bls12_381>> = (0..n).map(|_| BlsSigner::new(rng)).collect();
        let signers_pks_g2_1: Vec<_> = signers_1.iter().map(|s| s.bls_pk_g2).collect();
        let signers_pks_g1_1: Vec<_> = signers_1.iter().map(|s| s.bls_pk_g1).collect();

        let curr_ss = pvss::Params::<Bls12_381>::new(signers_pks_g2_0, t).unwrap();
        let next_ss = pvss::Params::<Bls12_381>::new(signers_pks_g2_1, t).unwrap();
        let bs_dkg = BsDkg::init(curr_ss, next_ss);

        let bs_transcript = bs_dkg.deal(signers_0[0].clone(), rng).unwrap();
        let verified_bs = bs_dkg.verify(bs_transcript, rng);
        let ss_1 = verified_bs.next_sharing;
        let ss_1_back = verified_bs.back_sharing;
        let tpk_1 = ThresholdVk::from_share(&ss_1.secret_sharing);

        // TWEAKS
        let h2_pred_1 = G2Affine::rand(rng); // TODO: hash_to_curve(C||1)
        let tweaks_1: Vec<_> = signers_1.iter()
            .map(|s| (s.sign_g2(h2_pred_1 - ss_1.secret_sharing.h2), s.bls_pk_g2))
            .collect();
        let tweaks_back_1: Vec<_> = signers_0.iter()
            .map(|s| (s.sign_g2(h2_pred_0 - ss_1_back.secret_sharing.h2), s.bls_pk_g2))
            .collect();

        let ss_mod_1 = ss_1.tweak_h2(&tweaks_1, h2_pred_1);
        let ss_back_mod_1 = ss_1_back.tweak_h2(&tweaks_back_1, h2_pred_0);

        let bgpk_delta = ss_mod_0.compute_delta(&ss_back_mod_1);

        let bgpk_mod_1: Vec<_> = ss_mod_1.tweaked_bgpks.iter().map(|x|x.unwrap()).collect();
        let tpk_mod_1 = Tpk {
            c: tpk_1.c,
            h2: h2_pred_1,
            g1: tpk_0.g1,
            g2: tpk_0.g2,
        };

        let sigs: Vec<_> = signers_1.iter().map(|s| s.sign_g1(message)).collect();
        let agg_sig_mod_1 = aggregate_sigs(sigs, signers_pks_g1_1, bgpk_mod_1, &config);
        tpk_mod_1.verify_sig(&agg_sig_mod_1, message, h2_pred_1);

        let tpk_1_pred: Tpk<Bls12_381> = Tpk {
            c: c_perm,
            h2: h2_pred_1,
            g1: tpk_0.g1,
            g2: tpk_0.g2,
        };
        let sig = Sig {
            asig: agg_sig_mod_1.asig,
            apk_g1: agg_sig_mod_1.apk_g1,
            apk_g2: agg_sig_mod_1.apk_g2,
            abgpk: (agg_sig_mod_1.abgpk - bgpk_delta).into_affine(),
        };
        tpk_1_pred.verify_sig(&sig, message, h2_pred_1);
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
