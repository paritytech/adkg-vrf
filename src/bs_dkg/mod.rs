use crate::bls::vanilla::BlsSigner;
use crate::dkg::transcript::{ContributionReceipt, Transcript};
use crate::pvss;
use crate::pvss::SecretSharing;
use crate::utils::BarycentricDomain;
use ark_ec::hashing::curve_maps::wb::{WBConfig, WBMap};
use ark_ec::hashing::map_to_curve_hasher::MapToCurve;
use ark_ec::pairing::Pairing;
use ark_ec::{AffineRepr, CurveGroup, VariableBaseMSM};
use ark_ff::Zero;
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

pub struct BsTranscript<C: Pairing> {
    //
    pub back_sharing: Transcript<C>,
    pub next_sharing: Transcript<C>,
}

pub struct BsDkg<C: Pairing> {
    pub curr: pvss::Params<C>,
    pub next: pvss::Params<C>,
}

impl<C: Pairing> BsDkg<C>
where
    <C::G2 as CurveGroup>::Config: WBConfig,
    WBMap<<C::G2 as CurveGroup>::Config>: MapToCurve<C::G2>,
{
    pub fn init(curr: pvss::Params<C>, next: pvss::Params<C>) -> Self {
        Self {
            curr,
            next,
        }
    }

    fn next(self, next: pvss::Params<C>) -> Self {
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
        let ss = self.next.deal_secrets(ssk, sh, rng)?;
        let receipt = ContributionReceipt::<C>::sign((ssk, ss.payload.c), (sh, ss.payload.h1), (dealer.sk, dealer.bls_pk_g1));
        let next_sharing = Transcript {
            agg_ss: ss,
            receipts: vec![(receipt, 1)],
        };

        let bs_sh = C::ScalarField::rand(rng);
        let bs_ss = self.curr.deal_secrets(ssk, bs_sh, rng)?;
        let bs_receipt = ContributionReceipt::<C>::sign((ssk, bs_ss.payload.c), (bs_sh, bs_ss.payload.h1), (dealer.sk, dealer.bls_pk_g1));
        let back_sharing = Transcript {
            agg_ss: bs_ss,
            receipts: vec![(bs_receipt, 1)],
        };

        Ok(BsTranscript { back_sharing, next_sharing })
    }

    fn verify<R: Rng>(&self, transcript: &BsTranscript<C>, rng: &mut R) {}

    pub fn tweak(ss: &SecretSharing<C>, tweaks: Vec<C::G2Affine>, h2_pred: C::G2Affine) -> BsKeys<C> {
        let tweaked_bgpks: Vec<C::G2> = ss.bgpk.iter()
            .zip(tweaks)
            .map(|(bgpk, tweak)| *bgpk + tweak)
            .collect();
        let tweaked_bgpks = C::G2::normalize_batch(&tweaked_bgpks);
        BsKeys {
            c: ss.c,
            bgpk: tweaked_bgpks,
            h2: h2_pred,
        }
    }

    pub fn compute_delta(&self, curr_keys: &BsKeys<C>, bs_next_keys: &BsKeys<C>) -> C::G2Affine {
        let bgpk_deltas: Vec<C::G2> = curr_keys.bgpk.iter()
            .zip(bs_next_keys.bgpk.iter())
            .map(|(curr, bs_next)| *bs_next - curr)
            .collect();
        let bgpk_deltas = C::G2::normalize_batch(&bgpk_deltas);
        let lis_at_zero = BarycentricDomain::of_size(self.curr.config.domain, self.curr.config.n)
            .lagrange_basis_at(C::ScalarField::zero());
        C::G2::msm(&bgpk_deltas, &lis_at_zero).unwrap()
            .into_affine()
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

