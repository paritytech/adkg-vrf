use crate::bls::vanilla::{verify_on_point, BlsSigInG2};
use crate::bs_dkg::sig_agg::EcSigAgg;
use crate::bs_dkg::{BsDkg, Committee};
use crate::dkg::transcript::Transcript;
use crate::hash_to_curve::PairingWithG2Map;
use crate::pvss::SecretSharing;
use crate::sig_agg::evaluate_at_0_in_g2;
use crate::pvss;
use ark_ec::pairing::Pairing;
use ark_ec::CurveGroup;
use ark_ff::Zero;
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
    tweaks: Vec<Option<C::G2Affine>>,

    /// The difference between the secrets of the current and the initial epochs.
    /// The secret key shared at epoch #K is `gsk_k = f_k(0).g2`.
    /// `self.gsk_delta = gsk_curr - gsk_0`, where `curr := self.sid`.
    /// Allows to verify signatures, produced using the shared secret key of the current epoch,
    /// against the public key `C = C_0 = f_0(0).g1` of the initial epoch.
    gsk_delta: C::G2,
}

impl<C: PairingWithG2Map> VerifiedSharingWithG1Keys<C> {
    pub fn from_silent_transcript(transcript: Transcript<C>, committee: &Committee<C>, sid: u64) -> Self {
        Self {
            sid,
            config: committee.params.config.clone(),
            ss: transcript.agg_ss.payload,
            signers_g1: committee.signers_g1.clone(),
            signers_g2: committee.params.signer_pks.clone(),
            h2_pred: BsDkg::<C>::h2_of(sid),
            tweaks: vec![],
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

impl<C: PairingWithG2Map> VerifiedSharingAndBack<C> {
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