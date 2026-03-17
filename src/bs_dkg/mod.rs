pub mod crypto;
pub mod sig_agg;

use crate::bls::vanilla::BlsSigner;
use crate::bs_dkg::sig_agg::EcSigAgg;
use crate::dkg::deal_and_sign;
use crate::dkg::transcript::{ContributionReceipt, Transcript};
use crate::hash_to_curve::PairingWithG2Map;
use crate::utils::BarycentricDomain;
use crate::{pvss, VerifiedSharing};
use ark_ec::pairing::Pairing;
use ark_ec::CurveGroup;
use ark_ec::VariableBaseMSM;
use ark_std::rand::Rng;
use ark_std::UniformRand;
use ark_std::Zero;


// TODO:
// 1. signers' pks in g1
// 2. signer's produce tweaks
// 3. tweaks are verified, aggregated and applied
// 4. new key material type
// 5. new signature type
// 6. aggregation for that
// basically 2 cryptosuites, ideally with a fall-back option

/// Back Sharing Distributed Key Generation Protocol.

/// The keys produced by the s-DKG protocol fix the committee of signers.
/// To enable advanced cryptography like encrypting to a future committee of signers
/// the bs-DKG keeps the threshold public key of an evolving set of signers predictable.

/// At each round a fresh secret is shared among a new set of signers,
/// and a related secret is (back-)shared to the current committee.
/// The new signers and the required threshold can be chosen arbitrarily.
/// After both the secrets shared are verified and aggregated,
/// the current and the next committee are each required to sign
/// (thus breaking the "silency" property of the original protocol).
/// These signatures together with the shared public (aggregation) keys allow to verify the
/// threshold signature of the new committee against the public key of the current committee.

/// The implementation assumes that the set of dealers for the next round is the current set of signers.

/// To make back-sharing useful, signers are assumed to publish their BLS public keys both in G1 and G2.
#[derive(Clone, Debug)]
pub struct Committee<C: Pairing> {
    /// Contains a list of signers' public keys in G2, specifies the threshold.
    pub params: pvss::Params<C>,
    /// Public keys of the same signers but in G1, in the same order, verified for consistency.
    pub signers_g1: Vec<C::G1Affine>,
}

pub struct BsTranscript<C: Pairing> {
    pub session_id: u64,
    pub back_sharing: Transcript<C>,
    pub next_sharing: Transcript<C>,
}

pub struct BsDkg<C: Pairing> {
    /// id of the next committee
    pub session_id: u64,
    /// Current set of signers (aka committee). Jointly know the secret of their epoch
    /// They share the new secret among the next set of signers,
    /// AND separately (via a different polynomial) "backshare" among themselves.
    /// Knowing the secrets of the 2 consecutive epochs, they
    /// collectively produce a key able to mutate threshold signatures produced by the different committees,
    /// that, in turn, keeps the public key persistent for a committee with evolving members.
    pub curr: Committee<C>,
    /// Next generation set of signers.
    pub next: Committee<C>,
}

impl<C: PairingWithG2Map> BsDkg<C> {
    pub fn start(next: Committee<C>) -> Self {
        Self {
            session_id: 0,
            curr: next.clone(),
            next,
        }
    }

    pub fn init(session_id: u64, curr: Committee<C>, next: Committee<C>) -> Self {
        Self {
            session_id,
            curr,
            next,
        }
    }

    pub fn next(self, next: Committee<C>) -> Self {
        Self {
            session_id: self.session_id + 1,
            curr: self.next,
            next,
        }
    }

    /// Predictable `h2` corresponding to the given round of the protocol.
    fn h2_of(session_id: u64) -> C::G2Affine {
        C::hash_to_g2(&session_id.to_be_bytes()).unwrap()
    }

    /// Predictable `h2` corresponding to the current round of the protocol.
    fn h2_curr(&self) -> C::G2Affine {
        Self::h2_of(self.session_id - 1)
    }

    /// Predictable `h2` corresponding to the next round of the protocol.
    fn h2_next(&self) -> C::G2Affine {
        Self::h2_of(self.session_id)
    }

    /// The same secret `ssk = f(0).g1` is
    /// 1. shared to the next committee with specified threshold `self.next.params.config.t`,
    /// 2. back-shared to the current committee with THEIR threshold `self.curr.params.config.t`.
    pub fn deal<R: Rng>(&self, dealer: BlsSigner<C>, rng: &mut R) -> Result<BsTranscript<C>, ()> {
        let ssk = C::ScalarField::rand(rng);

        let sh = C::ScalarField::rand(rng);
        let ss = self.next.params.deal_secrets(ssk, sh, rng)?;
        let receipt = ContributionReceipt::<C>::sign((ssk, ss.payload.c), (sh, ss.payload.h1), (dealer.sk, dealer.bls_pk_g1));
        let next_sharing = Transcript {
            agg_ss: ss,
            receipts: vec![(receipt, 1)],
        };

        let bs_sh = C::ScalarField::rand(rng);
        let bs_ss = self.curr.params.deal_secrets(ssk, bs_sh, rng)?;
        let bs_receipt = ContributionReceipt::<C>::sign((ssk, bs_ss.payload.c), (bs_sh, bs_ss.payload.h1), (dealer.sk, dealer.bls_pk_g1));
        let back_sharing = Transcript {
            agg_ss: bs_ss,
            receipts: vec![(bs_receipt, 1)],
        };

        Ok(BsTranscript { session_id: self.session_id, back_sharing, next_sharing })
    }

    pub fn deal_first<R: Rng>(&self, dealer: BlsSigner<C>, rng: &mut R) -> Result<Transcript<C>, ()> {
        deal_and_sign(&self.curr.params, rng, (dealer.sk, dealer.bls_pk_g1))
    }

    // TODO: this is a stub
    pub fn verify<R: Rng>(&self, bs_transcript: BsTranscript<C>, rng: &mut R) -> VerifiedSharingAndBack<C> {
        let BsTranscript {
            session_id,
            back_sharing,
            next_sharing,
        } = bs_transcript;
        assert_eq!(session_id, self.session_id);
        let back_sharing = VerifiedSharing {
            secret_sharing: back_sharing.agg_ss.payload,
            params: self.curr.params.clone(),
        };
        let back_sharing = VerifiedSharingWithG1Keys {
            verified_sharing: back_sharing,
            signers_g1: self.curr.signers_g1.clone(),
            h2_pred: self.h2_curr(),
        };
        let next_sharing = VerifiedSharing {
            secret_sharing: next_sharing.agg_ss.payload,
            params: self.next.params.clone(),
        };
        let next_sharing = VerifiedSharingWithG1Keys {
            verified_sharing: next_sharing,
            signers_g1: self.next.signers_g1.clone(),
            h2_pred: self.h2_next(),
        };
        VerifiedSharingAndBack {
            back_sharing,
            next_sharing,
        }
    }

    // TODO: this is a stub
    pub fn verify_first<R: Rng>(&self, s_transcript: Transcript<C>, rng: &mut R) -> VerifiedSharingWithG1Keys<C> {
        let verified_sharing = VerifiedSharing {
            secret_sharing: s_transcript.agg_ss.payload,
            params: self.curr.params.clone(),
        };
        VerifiedSharingWithG1Keys {
            verified_sharing,
            signers_g1: self.curr.signers_g1.clone(),
            h2_pred: self.h2_next(),
        }
    }
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

#[derive(Clone, Debug)]
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

    pub fn into_signature_aggregator(self) -> EcSigAgg<C> {
        EcSigAgg::new(self.secret_sharing.params.signer_pks, self.signers_g1, self.tweaked_bgpks, self.secret_sharing.params.config)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bls::threshold::ThresholdVk;
    use crate::bs_dkg::crypto::{EvolvingCommitteeAggSig, EvolvingCommitteePk, EvolvingCommitteeSig};
    use ark_bls12_381::Bls12_381;
    use ark_ec::AffineRepr;
    use ark_std::test_rng;


    fn committee<C: Pairing>(signers: &[BlsSigner<C>], t: usize) -> Committee<C> {
        let pks_g1: Vec<_> = signers.iter().map(|s| s.bls_pk_g1).collect();
        let pks_g2: Vec<_> = signers.iter().map(|s| s.bls_pk_g2).collect();
        Committee {
            params: pvss::Params::<C>::new(pks_g2, t).unwrap(),
            signers_g1: pks_g1,
        }
    }

    #[test]
    fn test_back_sharing() {
        let rng = &mut test_rng();

        // Protocol and the participants
        let (n, t) = (7, 5);

        let signers_0: Vec<BlsSigner<Bls12_381>> = (0..n).map(|_| BlsSigner::new(rng)).collect();
        let signers_1: Vec<BlsSigner<Bls12_381>> = (0..n).map(|_| BlsSigner::new(rng)).collect();
        let committee_0 = committee(&signers_0, t);
        let committee_1 = committee(&signers_1, t);
        let dealer = signers_0[0].clone();

        let bs_dkg = BsDkg::start(committee_0);

        // Deals secret shares to the epoch #1 committee (no-one to backshare to)
        let transcript = bs_dkg.deal_first(dealer.clone(), rng).unwrap();
        let ss_0 = bs_dkg.verify_first(transcript, rng);
        let h2_pred_0 = ss_0.h2_pred;
        let (fc_tpk_0, fc_sig_agg_0) = ss_0.verified_sharing.clone().into_keys();
        let ec_tpk = EvolvingCommitteePk::with_c(fc_tpk_0.c);

        // Tests a threshold signature at epoch 0
        let msg = BlsSigner::<Bls12_381>::hash_to_g1("message".as_bytes()).into_group();
        let sigs: Vec<_> = signers_0[..t].iter().map(|s| s.sign_g1(msg)).collect();
        let fc_asig_0 = fc_sig_agg_0.aggregate_wo_checking(sigs.clone());
        fc_tpk_0.verify_unoptimized(&fc_asig_0, msg);

        // Tweaks the key material (`bgpks` and `h2`)
        let tweak_msg_0 = ss_0.tweak_msg();
        let tweaks_0: Vec<_> = signers_0.iter()
            .map(|s| (s.sign_g2(tweak_msg_0), s.bls_pk_g2))
            .collect();
        let ss_tweaked_0 = ss_0.tweak(&tweaks_0);

        let ec_sig_agg_0 = ss_tweaked_0.clone().into_signature_aggregator();
        let ec_asig_0 = ec_sig_agg_0.aggregate(sigs);
        ec_tpk.verify_sig(&ec_asig_0, msg.into_affine(), h2_pred_0);

        let bs_dkg = bs_dkg.next(committee_1);
        let bs_transcript = bs_dkg.deal(dealer, rng).unwrap();
        let verified_bs = bs_dkg.verify(bs_transcript, rng);
        let ss_1 = verified_bs.next_sharing;
        let ss_1_back = verified_bs.back_sharing;
        let fc_tpk_1 = ThresholdVk::from_share(&ss_1.verified_sharing.secret_sharing);
        let ec_tpk_1 = EvolvingCommitteePk::with_c(fc_tpk_1.c);
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

        let ss_tweaked_1 = ss_1.tweak(&tweaks_1);
        let ss_back_tweaked_1 = ss_1_back.tweak(&tweaks_back_1);
        let bgpk_delta = ss_tweaked_0.compute_delta(&ss_back_tweaked_1);

        let sigs: Vec<_> = signers_1[n-t..].iter().map(|s| s.sign_g1(msg)).collect();
        let ec_sig_agg_1 = ss_tweaked_1.into_signature_aggregator();
        let ec_asig_1 = ec_sig_agg_1.aggregate(sigs);
        ec_tpk_1.verify_sig(&ec_asig_1, msg.into_affine(), h2_pred_1);

        let ec_asig_1 = EvolvingCommitteeAggSig(EvolvingCommitteeSig {
            sig: ec_asig_1.0.sig,
            pk_g1: ec_asig_1.0.pk_g1,
            pk_g2: ec_asig_1.0.pk_g2,
            bgpk: (ec_asig_1.0.bgpk - bgpk_delta).into_affine(),
        });
        ec_tpk.verify_sig(&ec_asig_1, msg.into_affine(), h2_pred_1);
    }
}