use crate::bls::vanilla::{verify_on_point, BlsSigInG2};
use crate::bs_dkg::sig_agg::EcSigAgg;
use crate::bs_dkg::{BsDkg, Committee};
use crate::dkg::transcript::Transcript;
use crate::hash_to_curve::PairingWithG2Map;
use crate::pvss;
use crate::pvss::SecretSharing;
use crate::sig_agg::evaluate_at_0_in_g2;
use ark_ec::pairing::Pairing;
use ark_ec::CurveGroup;
use ark_ff::Zero;
use ark_std::iterable::Iterable;
use hashbrown::HashMap;

/// Verified aggregated (related) secrets shared to `2` consequent committees of signers.
#[derive(Clone, Debug)]
pub struct VerifiedSharingAndBack<C: Pairing> {
    pub next_sharing: VerifiedSharingWithG1Keys<C>,
    pub back_sharing: VerifiedSharingWithG1Keys<C>,
}

#[derive(Clone, Debug)]
pub struct VerifiedSharingWithG1Keys<C: Pairing> {
    sid: u64,
    config: pvss::Config<C>,
    pub(crate) ss: SecretSharing<C>,
    signers_g1: Vec<C::G1Affine>,
    signers_g2: Vec<C::G2Affine>,
    h2_pred: C::G2Affine,
    /// `bgpks` adjusted to account the transition from `self.ss.h2` to `h2_pred`.
    /// `bgpk_j_tweaked = bgpk_j + tweak_j`
    /// `bgpk_j = f(w^j).g2 + sh.pk_j`
    /// `tweak_j = sk_j.(h2_pred - h2)`
    /// Thus `bgpk_j_tweaked = f(w^j).g2 + sk_j.h2_pred`.
    bgpks_tweaked: Vec<Option<C::G2Affine>>,

    /// The difference between the secrets of the current and the initial epochs.
    /// The secret key shared at epoch #K is `gsk_k = f_k(0).g2`.
    /// `self.gsk_delta = gsk_curr - gsk_0`, where `curr := self.sid`.
    /// Allows to verify signatures, produced using the shared secret key of the current epoch,
    /// against the public key `C = C_0 = f_0(0).g1` of the initial epoch.
    gsk_delta: C::G2,
}

impl<C: PairingWithG2Map> VerifiedSharingWithG1Keys<C> {
    pub fn from_silent_transcript(transcript: Transcript<C>, committee: &Committee<C>, sid: u64) -> Self {
        let config = committee.params.config.clone();
        Self {
            sid,
            config: config.clone(),
            ss: transcript.agg_ss.payload,
            signers_g1: committee.signers_g1.clone(),
            signers_g2: committee.params.signer_pks.clone(),
            h2_pred: BsDkg::<C>::h2_of(sid),
            bgpks_tweaked: vec![None; config.n],
            gsk_delta: C::G2::zero(),
        }
    }

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
        self.h2_pred - self.ss.h2
    }

    pub fn apply_tweaks(&mut self, tweak_sigs: Vec<BlsSigInG2<C>>) {
        let msg = self.h2_tweak().into_affine();
        // Signers are identified by their BLS public key in G1 sent together with the signature.
        let pk_to_j = self.pk_g1_to_j();
        // Remove tweaks that
        // 1. has been already applied,
        // 2. `pk` from the signature doesn't belong to the signer set.
        let useful_sigs: Vec<(usize, BlsSigInG2<C>)> = tweak_sigs.into_iter().filter_map(|sig| {
            pk_to_j.get(&sig.pk) // signature's pk recognized
                .and_then(|&j| self.bgpks_tweaked[j].is_none() // `bgpk` isn't tweaked yet
                    .then(|| (j, sig)))
        }).collect();
        // TODO: BLS batch verification
        let valid_sigs: Vec<(usize, BlsSigInG2<C>)> = useful_sigs.into_iter()
            .filter(|(_, sig)| verify_on_point::<C>(sig.sig, msg, sig.pk, self.config.g1))
            .collect();
        // `bgpk_j_tweaked = bgpk_j + tweak_j`
        let (js, tweaked_bgpks): (Vec<usize>, Vec<C::G2>) = valid_sigs.into_iter()
            .map(|(j, sig)| (j, self.ss.bgpk[j] + sig.sig))
            .unzip();
        let tweaked_bgpks = C::G2::normalize_batch(&tweaked_bgpks);
        js.into_iter().zip(tweaked_bgpks)
            .for_each(|(j, bgpk_j)| self.bgpks_tweaked[j] = Some(bgpk_j));
    }

    fn pk_g1_to_j(&self) -> HashMap<C::G1Affine, usize> {
        self.signers_g1.iter()
            .copied()
            .enumerate()
            .map(|(j, pk_g1)| (pk_g1, j))
            .collect()
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
    /// Thus, to compute the `delta`, there should be not less than `t0 = self.config.t`
    /// signers who have their `bgpk`s tweaked in bot the current share and the back share.
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
    fn compute_next_gsk_delta(&self, bs: &Self) -> C::G2 {
        let epoch_delta = self.compute_epoch_gsk_delta(bs);
        let delta_curr = self.gsk_delta;
        let delta_next = delta_curr + epoch_delta;
        delta_next
    }

    // Returns `Some((curr.tweaked_bgpk_j, bs.tweaked_bgpk_j))` at the `j`-th position, or `None`.
    fn get_tweaked_bgpk_pairs(&self, bs: &Self) -> impl Iterator<Item=Option<(C::G2Affine, C::G2Affine)>> {
        self.bgpks_tweaked.clone().into_iter()
            .zip(bs.bgpks_tweaked.clone().into_iter())
            .map(|(curr_opt, back_opt)|
                     curr_opt.zip(back_opt) // both `Options` should be `Some`
            )
    }

    pub fn tweaked_signers_count(&self) -> usize {
        self.bgpks_tweaked.iter().flatten().count()
    }

    fn can_interpolate_delta(&self, bs: &Self) -> bool {
        self.get_tweaked_bgpk_pairs(bs).flatten().count() >= self.config.t
    }

    fn can_interpolate_anything(&self) -> bool {
        self.tweaked_signers_count() >= self.config.t
    }

    pub fn into_combiner(self) -> Option<EcSigAgg<C>> {
        self.can_interpolate_anything().then(||
            EcSigAgg::new(self.sid,
                          self.signers_g2,
                          self.signers_g1,
                          self.bgpks_tweaked,
                          self.gsk_delta,
                          self.config,
            ))
    }
}

impl<C: PairingWithG2Map> VerifiedSharingAndBack<C> {

    /// Updates the `gsk_delta` of the next committee.
    /// It is required to verify threshold proofs produced by the next committee with the
    /// permanent public key (the public key of the committee at epoch #0).
    /// `next.gsk_delta = curr.gsk_delta + delta(curr, back)`.
    ///
    /// Consumes `self.back_sharing` as it was only need to compute this delta.
    fn sharing_with_delta(self, curr: &VerifiedSharingWithG1Keys<C>) -> VerifiedSharingWithG1Keys<C> {
        let mut next = self.next_sharing;
        next.gsk_delta = curr.compute_next_gsk_delta(&self.back_sharing);
        next
    }

    pub fn tweak_bs_and_try_compute_delta(mut self, tweak_sigs: Vec<BlsSigInG2<C>>, curr: &VerifiedSharingWithG1Keys<C>) -> Option<VerifiedSharingWithG1Keys<C>> {
        self.back_sharing.apply_tweaks(tweak_sigs);
        curr.can_interpolate_delta(&self.back_sharing).then(|| self.sharing_with_delta(curr))
    }
}