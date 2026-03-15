use crate::bls::vanilla::BlsSigner;
use crate::dkg::deal_and_sign;
use crate::dkg::transcript::{ContributionReceipt, Transcript};
use crate::{pvss, VerifiedSharing, VerifiedSharingAndBack, VerifiedSharingWithG1Keys};
use ark_ec::hashing::curve_maps::wb::{WBConfig, WBMap};
use ark_ec::hashing::map_to_curve_hasher::MapToCurve;
use ark_ec::pairing::Pairing;
use ark_ec::{AffineRepr, CurveGroup};
use ark_std::rand::Rng;
use ark_std::UniformRand;

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
pub struct Committee<C: Pairing> {
    /// Contains a list of signers' public keys in G2, specifies the threshold.
    pub params: pvss::Params<C>,
    /// Public keys of the same signers but in G1, in the same order, verified for consistency.
    pub signers_g1: Vec<C::G1Affine>,
}

pub struct BsTranscript<C: Pairing> {
    pub back_sharing: Transcript<C>,
    pub next_sharing: Transcript<C>,
}

pub struct BsDkg<C: Pairing> {
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

impl<C: Pairing> BsDkg<C>
where
    <C::G2 as CurveGroup>::Config: WBConfig,
    WBMap<<C::G2 as CurveGroup>::Config>: MapToCurve<C::G2>,
{
    pub fn init(curr: Committee<C>, next: Committee<C>) -> Self {
        Self {
            curr,
            next,
        }
    }

    fn next(self, next: Committee<C>) -> Self {
        Self {
            curr: self.next,
            next,
        }
    }

    fn h2() -> C::G2Affine {
        C::G2Affine::generator() // TODO
    }

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

        Ok(BsTranscript { back_sharing, next_sharing })
    }

    pub fn deal_first<R: Rng>(&self, dealer: BlsSigner<C>, rng: &mut R) -> Result<Transcript<C>, ()> {
        deal_and_sign(&self.curr.params, rng, (dealer.sk, dealer.bls_pk_g1))
    }

    // TODO: this is a stub
    pub fn verify<R: Rng>(&self, bs_transcript: BsTranscript<C>, rng: &mut R) -> VerifiedSharingAndBack<C> {
        let BsTranscript {
            back_sharing,
            next_sharing,
        } = bs_transcript;
        let back_sharing = VerifiedSharing {
            secret_sharing: back_sharing.agg_ss.payload,
            params: self.curr.params.clone(),
        };
        let back_sharing = VerifiedSharingWithG1Keys {
            verified_sharing: back_sharing,
            signers_g1: self.curr.signers_g1.clone(),
            h2_pred: C::G2Affine::rand(rng),
        };
        let next_sharing = VerifiedSharing {
            secret_sharing: next_sharing.agg_ss.payload,
            params: self.next.params.clone(),
        };
        let next_sharing = VerifiedSharingWithG1Keys {
            verified_sharing: next_sharing,
            signers_g1: self.curr.signers_g1.clone(),
            h2_pred: C::G2Affine::rand(rng),
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
            h2_pred: C::G2Affine::rand(rng), // TODO: make predictable
        }
    }
}

pub struct BsKeys<C: Pairing> {
    /// The public key corresponding to the shared secret key.
    /// `c = f(0).g1`
    pub c: C::G1Affine,
    /// Shares of the secret, encrypted to the signers.
    /// `(bgpk_j = f(w^j).g2 + sh.pk_j, j = 0,...,n-1)`
    pub bgpk: Vec<C::G2Affine>,
    /// `h2 = sh.g2`
    pub h2: C::G2Affine,
}

