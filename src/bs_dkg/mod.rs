pub mod crypto;
pub mod sig_agg;

use crate::bls::vanilla::{verify_on_point, BlsSigInG2, BlsSigner};
use crate::bs_dkg::sig_agg::EcSigAgg;
use crate::dkg::transcript::Transcript;
use crate::dkg::{deal_and_sign, deal_and_sign_ssk};
use crate::hash_to_curve::PairingWithG2Map;
use crate::pvss::SecretSharing;
use crate::sig_agg::prepare;
use crate::{pvss, VerifiedSharing};
use ark_ec::pairing::Pairing;
use ark_ec::VariableBaseMSM;
use ark_ec::{AffineRepr, CurveGroup};
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
            abgpk_delta: C::G2Affine::zero(),
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
            abgpk_delta: C::G2Affine::zero(),
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
            abgpk_delta: C::G2Affine::zero(),
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
    abgpk_delta: C::G2Affine,
}

impl<C: Pairing> VerifiedSharingWithG1Keys<C> {
    pub fn tweak_msg(&self) -> C::G2 {
        self.h2_pred - self.verified_sharing.secret_sharing.h2
    }

    // pub fn tweak_bgpks(&self, sigs: &[(C::G2Affine, C::G2Affine)]) -> Vec<Option<C::G2Affine>> {
    //     let tweak_msg = self.tweak_msg().into_affine();
    //     let tweaked_bgpks = self.verified_sharing.aggregation_key().aggregate_tweaks(tweak_msg, sigs);
    //     tweaked_bgpks
    // }

    // pub fn tweak(self, sigs: &[(C::G2Affine, C::G2Affine)]) -> TweakedSharing<C> {
    //     // TODO: verify
    //     let tweaked_bgpks = self.tweak_bgpks(sigs);
    //     TweakedSharing {
    //         sid: self.sid,
    //         secret_sharing: self.verified_sharing,
    //         signers_g1: self.signers_g1,
    //         tweaked_bgpks,
    //         h2_pred: self.h2_pred,
    //     }
    // }

    pub fn tweak_self(&mut self, sigs: Vec<BlsSigInG2<C>>) {
        let tweaks = self.prepare_tweaks(sigs);
        self.tweaks = tweaks; //TODO: merge?
    }

    pub fn prepare_tweaks(&self, sigs: Vec<BlsSigInG2<C>>) -> Vec<Option<C::G2Affine>> {
        let mut tweaks = vec![None; self.config.n];
        let msg = self.tweak_msg().into_affine();
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

    /// Computes (f0(0) - f1_back(0)).g2
    /// `deg(f0) = deg(f1_back) = t0 - 1`.
    pub fn compute_delta(&self, bs: &Self) -> C::G2Affine {
        let bgpk_deltas: Vec<Option<C::G2>> = self.tweaked_bpks().into_iter()
            .zip(bs.tweaked_bpks().into_iter())
            .map(|(curr, bs)| curr.zip(bs).map(|(curr, bs)| bs - curr))
            .collect();
        let (lis_at_zero, bgpk_deltas) = prepare(bgpk_deltas, &self.config);
        let bgpk_deltas = C::G2::normalize_batch(&bgpk_deltas);
        C::G2::msm(&bgpk_deltas, &lis_at_zero).unwrap().into_affine()
    }

    pub fn into_combiner(self) -> EcSigAgg<C> {
        let bgpks = self.tweaks.iter().enumerate().map(|(j, sig)| sig.map(|sig| (sig + self.ss.bgpk[j]).into_affine())).collect();
        EcSigAgg::new(self.sid, self.signers_g2, self.signers_g1, bgpks, self.config)
    }
}


/// Verified aggregated (related) secrets shared to `2` consequent committees of signers.
#[derive(Clone, Debug)]
pub struct VerifiedSharingAndBack<C: Pairing> {
    back_sharing: VerifiedSharingWithG1Keys<C>,
    next_sharing: VerifiedSharingWithG1Keys<C>,
}

// #[derive(Clone, Debug)]
// pub struct TweakedSharing<C: Pairing> {
//     sid: u64,
//     /// Secret shared to signers with keys in G2.
//     secret_sharing: VerifiedSharing<C>,
//     /// Signers' keys in G1 in the same order as `self.secret_sharing.params.signer_pks`
//     signers_g1: Vec<C::G1Affine>,
//     tweaked_bgpks: Vec<Option<C::G2Affine>>,
//     h2_pred: C::G2Affine,
// }
//
// impl<C: Pairing> TweakedSharing<C> {
//     /// Computes the threshold public key delta. This is (f0(0) - f1_back(0)).g1
//     /// `deg(f0) = deg(f1_back) = t0 - 1`.
//     /// Thus both the current sharing and the back-shared one should have `t0` tweaked bgpks at the same positions.
//     /// TODO: check for that
//     pub fn compute_delta(&self, bs_next: &Self) -> C::G2Affine {
//         let bgpk_deltas: Vec<C::G2> = self.tweaked_bgpks.iter()
//             .zip(bs_next.tweaked_bgpks.iter())
//             .map(|(curr, bs_next)| bs_next.unwrap() - curr.unwrap())
//             .collect();
//         let bgpk_deltas = C::G2::normalize_batch(&bgpk_deltas);
//         let bs_next_config = &bs_next.secret_sharing.params.config;
//         let lis_at_zero = BarycentricDomain::of_size(bs_next_config.domain, bs_next_config.n)
//             .lagrange_basis_at(C::ScalarField::zero());
//         C::G2::msm(&bgpk_deltas, &lis_at_zero).unwrap()
//             .into_affine()
//     }
//
//     pub fn into_signature_aggregator(self) -> EcSigAgg<C> {
//         EcSigAgg::new(self.sid, self.secret_sharing.params.signer_pks, self.signers_g1, self.tweaked_bgpks, self.secret_sharing.params.config)
//     }
// }

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
        let tweak_msg_0 = ss_0.tweak_msg();
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
        let verified_bs = bs_dkg.verify(bs_transcript, rng);
        let mut ss_1 = verified_bs.next_sharing;
        let c1 = ss_1.verified_sharing.secret_sharing.c;
        let mut ss_1_back = verified_bs.back_sharing;

        // TWEAKS
        let tweak_msg_1 = ss_1.tweak_msg();
        let tweaks_1: Vec<_> = signers_1[n - t..].iter()
            .map(|s| s.sign_g2_point(tweak_msg_1))
            .collect();
        let tweak_msg_back_1 = ss_1_back.tweak_msg();
        let tweaks_back_1: Vec<_> = signers_0[..t].iter()
            .map(|s| s.sign_g2_point(tweak_msg_back_1))
            .collect();
        // let ss_tweaked_1 = ss_1.tweak(&tweaks_1);
        ss_1.tweak_self(tweaks_1);
        // let ss_back_tweaked_1 = ss_1_back.tweak(&tweaks_back_1);
        ss_1_back.tweak_self(tweaks_back_1);
        let bgpk_delta = ss_0.compute_delta(&ss_1_back);
        let ec_pk_1 = EvolvingCommitteePk::with_c(c1);

        let sig_agg_1 = ss_1.into_combiner();
        let sigs_1: Vec<_> = signers_1[n - t..].iter().map(|s| s.sign_bytes_in_g1(b"msg1")).collect();
        let mut asig_1 = sig_agg_1.aggregate(sigs_1);
        ec_pk_1.verify(&asig_1, b"msg1");

        asig_1.asig.bgpk = (asig_1.asig.bgpk - bgpk_delta).into_affine();
        ec_pk.verify(&asig_1, b"msg1");
    }
}