pub mod crypto;
pub mod sig_agg;
mod sharing;

use crate::bls::vanilla::BlsSigner;
use crate::dkg::transcript::Transcript;
use crate::dkg::{deal_and_sign, deal_and_sign_ssk};
use crate::hash_to_curve::PairingWithG2Map;
use crate::pvss;
use ark_ec::pairing::Pairing;
use ark_std::rand::Rng;
use ark_std::UniformRand;
use sharing::{VerifiedSharingAndBack, VerifiedSharingWithG1Keys};
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
    pub sid: u64,
    pub next_transcript: Transcript<C>,
    pub back_transcript: Transcript<C>,
}

pub struct BsDkg<C: Pairing> {
    /// ordinal of the next committee
    pub sid: u64,
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
            sid: 0,
            curr: next.clone(),
            next,
        }
    }

    pub fn init(sid: u64, curr: Committee<C>, next: Committee<C>) -> Self {
        Self {
            sid,
            curr,
            next,
        }
    }

    pub fn next(self, next: Committee<C>) -> Self {
        Self {
            sid: self.sid + 1,
            curr: self.next,
            next,
        }
    }

    /// Predictable `h2` corresponding to the given round of the protocol.
    fn h2_of(sid: u64) -> C::G2Affine {
        C::hash_to_g2(&sid.to_be_bytes()).unwrap()
    }

    /// Predictable `h2` corresponding to the current round of the protocol.
    fn h2_curr(&self) -> C::G2Affine {
        Self::h2_of(self.sid - 1)
    }

    /// Predictable `h2` corresponding to the next round of the protocol.
    fn h2_next(&self) -> C::G2Affine {
        Self::h2_of(self.sid)
    }

    /// The same secret `ssk = f(0).g1` is
    /// 1. shared to the next committee with specified threshold `self.next.params.config.t`,
    /// 2. back-shared to the current committee with THEIR threshold `self.curr.params.config.t`.
    pub fn deal<R: Rng>(&self, dealer: BlsSigner<C>, rng: &mut R) -> Result<BsTranscript<C>, ()> {
        let ssk = C::ScalarField::rand(rng);
        let next_transcript = deal_and_sign_ssk(ssk, &self.next.params, rng, dealer.as_tuple()).unwrap();
        let back_transcript = deal_and_sign_ssk(ssk, &self.curr.params, rng, dealer.as_tuple()).unwrap();
        Ok(BsTranscript { sid: self.sid, next_transcript, back_transcript })
    }

    pub fn deal_first<R: Rng>(&self, dealer: BlsSigner<C>, rng: &mut R) -> Result<Transcript<C>, ()> {
        deal_and_sign(&self.curr.params, rng, (dealer.sk, dealer.pk_g1))
    }

    // TODO: this is a stub
    pub fn verify<R: Rng>(&self, bs_transcript: BsTranscript<C>, rng: &mut R) -> VerifiedSharingAndBack<C> {
        let BsTranscript {
            sid,
            next_transcript,
            back_transcript,
        } = bs_transcript;
        assert_eq!(sid, self.sid);
        let next_sharing = VerifiedSharingWithG1Keys::from_silent_transcript(next_transcript, &self.next, self.sid);
        let back_sharing = VerifiedSharingWithG1Keys::from_silent_transcript(back_transcript, &self.curr, self.sid - 1);
        VerifiedSharingAndBack {
            back_sharing,
            next_sharing,
        }
    }

    // TODO: this is a stub
    pub fn verify_first<R: Rng>(&self, transcript: Transcript<C>, rng: &mut R) -> VerifiedSharingWithG1Keys<C> {
        VerifiedSharingWithG1Keys::from_silent_transcript(transcript, &self.curr, self.sid)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::bs_dkg::crypto::EvolvingCommitteePk;
    use ark_bls12_381::Bls12_381;
    use ark_std::test_rng;

    fn committee<C: Pairing>(signers: &[BlsSigner<C>], t: usize) -> Committee<C> {
        let (pks_g1, pks_g2): (Vec<_>, Vec<_>) = signers.iter().map(|s| (s.pk_g1, s.pk_g2)).unzip();
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

        // EPOCH #0

        // Deals secret shares to the epoch #0 committee (noone to backshare to)
        let bs_dkg = BsDkg::start(committee_0);
        let transcript = bs_dkg.deal_first(dealer.clone(), rng).unwrap();
        let mut ss_0 = bs_dkg.verify_first(transcript, rng);
        let c0 = ss_0.ss.c;
        let ec_pk = EvolvingCommitteePk::with_c(c0);

        // Tweaks the key material (`bgpks` and `h2`)
        let tweak_msg_0 = ss_0.h2_tweak();
        let tweaks_0: Vec<_> = signers_0[..t].iter().map(|s| s.sign_g2_point(tweak_msg_0)).collect();
        // let ss_tweaked_0 = ss_0.tweak(&tweaks_0);
        ss_0.tweak_self(tweaks_0);

        // Tests a threshold signature at epoch 0
        let sig_agg_0 = ss_0.clone().into_combiner();
        let sigs_0: Vec<_> = signers_0[..t].iter().map(|s| s.sign_bytes_in_g1(b"msg0")).collect();
        let asig_0 = sig_agg_0.aggregate(sigs_0);
        ec_pk.verify(&asig_0, b"msg0");

        // EPOCH #1
        let bs_dkg = bs_dkg.next(committee_1);
        let bs_transcript = bs_dkg.deal(dealer, rng).unwrap();
        let mut verified_bs = bs_dkg.verify(bs_transcript, rng);
        let ss_1 = &mut verified_bs.next_sharing;
        let c1 = ss_1.ss.c;
        let ss_1_back = &mut verified_bs.back_sharing;

        // TWEAKS
        let tweak_msg_1 = ss_1.h2_tweak();
        let tweaks_1: Vec<_> = signers_1[n - t..].iter()
            .map(|s| s.sign_g2_point(tweak_msg_1))
            .collect();
        let tweak_msg_back_1 = ss_1_back.h2_tweak();
        let tweaks_back_1: Vec<_> = signers_0[..t].iter()
            .map(|s| s.sign_g2_point(tweak_msg_back_1))
            .collect();
        // let ss_tweaked_1 = ss_1.tweak(&tweaks_1);
        ss_1.tweak_self(tweaks_1);
        // let ss_back_tweaked_1 = ss_1_back.tweak(&tweaks_back_1);
        ss_1_back.tweak_self(tweaks_back_1);
        // let bgpk_delta = ss_0.compute_delta(&ss_1_back);
        let ss_1 = verified_bs.sharing_with_delta(&ss_0);

        let sig_agg_1 = ss_1.into_combiner();
        let sigs_1: Vec<_> = signers_1[n - t..].iter().map(|s| s.sign_bytes_in_g1(b"msg1")).collect();
        let asig_1 = sig_agg_1.aggregate(sigs_1);
        ec_pk.verify(&asig_1, b"msg1");
    }
}