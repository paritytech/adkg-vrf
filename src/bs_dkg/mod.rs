pub mod crypto;
pub mod sig_agg;

use crate::bls::vanilla::{verify_on_point, BlsSigInG2, BlsSigner};
use crate::bs_dkg::sig_agg::EcSigAgg;
use crate::dkg::transcript::Transcript;
use crate::dkg::{deal_and_sign, deal_and_sign_ssk};
use crate::hash_to_curve::PairingWithG2Map;
use crate::pvss::SecretSharing;
use crate::sig_agg::evaluate_at_0_in_g2;
use crate::{pvss, VerifiedSharing};
use ark_ec::pairing::Pairing;
use ark_ec::CurveGroup;
use ark_ff::Zero;
use ark_std::rand::Rng;
use ark_std::UniformRand;
use hashbrown::HashMap;
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
        let next_sharing = deal_and_sign_ssk(ssk, &self.next.params, rng, dealer.as_tuple()).unwrap();
        let back_sharing = deal_and_sign_ssk(ssk, &self.curr.params, rng, dealer.as_tuple()).unwrap();
        Ok(BsTranscript { session_id: self.session_id, back_sharing, next_sharing })
    }

    pub fn deal_first<R: Rng>(&self, dealer: BlsSigner<C>, rng: &mut R) -> Result<Transcript<C>, ()> {
        deal_and_sign(&self.curr.params, rng, (dealer.sk, dealer.pk_g1))
    }

    // TODO: this is a stub
    pub fn verify<R: Rng>(&self, bs_transcript: BsTranscript<C>, rng: &mut R) -> VerifiedSharingAndBack<C> {
        let BsTranscript {
            session_id,
            back_sharing,
            next_sharing,
        } = bs_transcript;
        assert_eq!(session_id, self.session_id);
        let _back_sharing = VerifiedSharing {
            secret_sharing: back_sharing.agg_ss.payload.clone(),
            params: self.curr.params.clone(),
        };
        let back_sharing = VerifiedSharingWithG1Keys {
            sid: session_id,
            verified_sharing: _back_sharing,
            config: self.curr.params.config.clone(),
            ss: back_sharing.agg_ss.payload,
            signers_g1: self.curr.signers_g1.clone(),
            signers_g2: self.curr.params.signer_pks.clone(),
            h2_pred: self.h2_curr(),
            tweaks: vec![],
            gsk_delta: C::G2::zero(),
        };
        let _next_sharing = VerifiedSharing {
            secret_sharing: next_sharing.agg_ss.payload.clone(),
            params: self.next.params.clone(),
        };
        let next_sharing = VerifiedSharingWithG1Keys {
            sid: session_id,
            config: self.next.params.config.clone(),
            ss: next_sharing.agg_ss.payload,
            verified_sharing: _next_sharing,
            signers_g1: self.next.signers_g1.clone(),
            signers_g2: self.next.params.signer_pks.clone(),
            h2_pred: self.h2_next(),
            tweaks: vec![],
            gsk_delta: C::G2::zero(),
        };

        VerifiedSharingAndBack {
            back_sharing,
            next_sharing,
        }
    }

    // TODO: this is a stub
    pub fn verify_first<R: Rng>(&self, s_transcript: Transcript<C>, rng: &mut R) -> VerifiedSharingWithG1Keys<C> {
        let verified_sharing = VerifiedSharing {
            secret_sharing: s_transcript.agg_ss.payload.clone(),
            params: self.curr.params.clone(),
        };
        VerifiedSharingWithG1Keys {
            sid: self.session_id,
            verified_sharing,
            config: self.curr.params.config.clone(),
            ss: s_transcript.agg_ss.payload,
            signers_g1: self.curr.signers_g1.clone(),
            signers_g2: self.curr.params.signer_pks.clone(),
            h2_pred: self.h2_next(),
            tweaks: vec![],
            gsk_delta: C::G2::zero(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct VerifiedSharingWithG1Keys<C: Pairing> {
    sid: u64,
    verified_sharing: VerifiedSharing<C>,
    config: pvss::Config<C>,
    ss: SecretSharing<C>,
    signers_g1: Vec<C::G1Affine>,
    signers_g2: Vec<C::G2Affine>,
    h2_pred: C::G2Affine,
    tweaks: Vec<Option<C::G2Affine>>,

    /// The difference between the secrets of the current and the initial epochs.
    /// The secret key shared at epoch #K is `gsk_k = f_k(0).g2`.
    /// `self.gsk_delta = gsk_curr - gsk_0`, where `curr := self.sid`.
    /// Allows to verify signatures, produced using the shared secret key of the current epoch,
    /// against the public key `C = C_0 = f_0(0).g1` of the initial epoch.
    gsk_delta: C::G2,
}

impl<C: Pairing> VerifiedSharingWithG1Keys<C> {

    /// The final sharing contains `(h1, h2)`, such that `h1 = sh.h1, h2 = sh.g2` for some `sh`,
    /// that is used to encrypt the `j`-the signer's share of the secret `f(w_j).g2` with ElGamal as
    /// `bgpk_j = f(w^j).g2 + sh.pk_j = f(w^j).g2 + sk_j.h2` for `sk_j.g2 = pk_j.`
    ///
    /// We want to replace this random `h2 = sh.g2` with a predictable `h2_perm = hash_to_g2(sid)`.
    /// For that each signer, is required to publish a BLS signature (in G2) on `h2_tweak = h2_pred - h2`.
    /// `tweak_j = sk_j.h2_tweak = sk_j.h2_pred - sk_j.h2`.
    ///
    /// That allows to adjust `bgpk_j` to `bgpk_j_tweaked = bgpk_j + tweak_j = f(w^j).g2 + sk_j.h2_pred` accordingly.
    ///
    /// If it doesn't happen signer `j`'s signatures can't be aggregated in the evolving committee scheme.
    ///
    pub fn h2_tweak(&self) -> C::G2 {
        self.h2_pred - self.verified_sharing.secret_sharing.h2
    }

    pub fn tweak_self(&mut self, sigs: Vec<BlsSigInG2<C>>) {
        let tweaks = self.prepare_tweaks(sigs);
        self.tweaks = tweaks; //TODO: merge?
    }

    pub fn prepare_tweaks(&self, sigs: Vec<BlsSigInG2<C>>) -> Vec<Option<C::G2Affine>> {
        let mut tweaks = vec![None; self.config.n];
        let msg = self.h2_tweak().into_affine();
        let pk_to_j = self.pk_g1_to_j();
        sigs.iter().for_each(|sig| {
            if let Some(&j) = pk_to_j.get(&sig.pk) {
                if verify_on_point::<C>(sig.sig, msg, sig.pk, self.config.g1) {
                    tweaks[j] = Some(sig.sig)
                }
            }
        });
        tweaks
    }

    pub fn pk_g1_to_j(&self) -> HashMap<C::G1Affine, usize> {
        self.signers_g1.iter()
            .copied()
            .enumerate()
            .map(|(j, pk_g1)| (pk_g1, j))
            .collect()
    }

    pub fn tweaked_bpks(&self) -> Vec<Option<C::G2>> {
        self.tweaks.iter().enumerate().map(|(j, sig)| sig.map(|sig| sig + self.ss.bgpk[j])).collect()
    }

    /// Computes the difference between the secrets `gsk = ssk.g2 = f(0).g2` of 2 subsequent committees
    ///  `delta = f_next(0).g2 - f_curr(0).g2`.
    ///
    /// Let `f0 := f_curr` be the current secret polynomial of degree `t0 - 1`.
    /// The degree of `f1 := f_next` is not relevant.
    ///
    /// After applying the tweaks shares of the current committee are
    /// `bgpk_j_tweaked = f0(w^j).g2 + sk_j.h2_pred_0, j=1,...,n0`
    ///
    /// The next epoch secret `f1(0).g2` is BACK-shared to the same (current) committee
    /// using a degree `t0 - 1` polynomial `f_back := f'1` such that `f'1(0) = f1(0)`.
    /// The back shares are tweaked to have the same `h2_pred_0`:
    /// `bgpk_back_j_tweaked = f'1(w^j).g2 + sk_j.h2_pred_0, j=1,...,n0`.
    ///
    /// Then `delta_j = bgpk_back_j_tweaked - bgpk_j_tweaked = (f'1(w^j) - f0(w^j)).g2`.
    /// The polynomial `f'1 - f0` has degree `t0 - 1`,
    /// so `(f'1 - f0)(0).g2` can be interpolated from `t0` such `delta_j`s.
    /// Finally, `f'1(0) - f0(0) = f1(0) - f0(0)`.
    ///
    /// So the requirement for this method to succeed is to have not less than `t0 = self.config.t`
    /// pairs of corresponding `bgpk`s tweaked in the current share and the BACK share.
    ///
    fn compute_epoch_gsk_delta(&self, bs: &Self) -> C::G2 {
        // `(f_back - f_curr)(w^j)` for some `j`s
        let f_deltas: Vec<Option<C::G2>> = self.get_tweaked_bgpk_pairs(bs)
            .map(|opt| opt.map(|(curr, back)| back - curr))
            .collect();
        let delta_f_at_0 = evaluate_at_0_in_g2(f_deltas, &self.config).unwrap(); //TODO
        delta_f_at_0
    }

    /// Computes the difference `delta_next = gsk_next - gsk_0` between the secrets of the next and the initial epochs
    ///
    /// `delta_curr := self.gsk_delta = gsk_curr - gsk_0`
    /// `epoch_delta = gsk_next - gsk_curr`
    /// `delta_next = gsk_next - gsk_0 = self.gsk_delta + epoch_delta`
    pub(crate) fn compute_next_gsk_delta(&self, bs: &Self) -> C::G2 {
        let epoch_delta = self.compute_epoch_gsk_delta(bs);
        let delta_curr = self.gsk_delta;
        let delta_next = delta_curr + epoch_delta;
        delta_next
    }

    // Returns `Some((curr.tweaked_bgpk_j, bs.tweaked_bgpk_j))` at the `j`-th position, or `None`.
    fn get_tweaked_bgpk_pairs(&self, bs: &Self) -> impl Iterator<Item=Option<(C::G2, C::G2)>> {
        self.tweaked_bpks().into_iter()
            .zip(bs.tweaked_bpks().into_iter())
            .map(|(curr_opt, back_opt)|
                     curr_opt.zip(back_opt) // both `Options` should be `Some`
            )
    }

    fn can_interpolate_delta(&self, bs: &Self) -> bool {
        self.get_tweaked_bgpk_pairs(bs).flatten().count() >= self.config.t
    }

    pub fn into_combiner(self) -> EcSigAgg<C> {
        let bgpks = self.tweaks.iter().enumerate().map(|(j, sig)| sig.map(|sig| (sig + self.ss.bgpk[j]).into_affine())).collect();
        EcSigAgg::new(self.sid, self.signers_g2, self.signers_g1, bgpks, self.gsk_delta, self.config)
    }
}


/// Verified aggregated (related) secrets shared to `2` consequent committees of signers.
#[derive(Clone, Debug)]
pub struct VerifiedSharingAndBack<C: Pairing> {
    back_sharing: VerifiedSharingWithG1Keys<C>,
    next_sharing: VerifiedSharingWithG1Keys<C>,
}

impl<C: Pairing> VerifiedSharingAndBack<C> {

    /// Updates the `gsk_delta` of the next committee.
    /// It is required to verify threshold proofs produced by the next committee with the
    /// permanent public key (the public key of the committee at epoch #0).
    /// `next.gsk_delta = curr.gsk_delta + delta(curr, back)`.
    ///
    /// Consumes `self.back_sharing` as it was only need to compute this delta.
    pub fn sharing_with_delta(self, curr: &VerifiedSharingWithG1Keys<C>) -> VerifiedSharingWithG1Keys<C> {
        let mut next = self.next_sharing;
        next.gsk_delta = curr.compute_next_gsk_delta(&self.back_sharing);
        next
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
        let ec_pk = EvolvingCommitteePk::with_c(ss_0.verified_sharing.secret_sharing.c);

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
        let c1 = ss_1.verified_sharing.secret_sharing.c;
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