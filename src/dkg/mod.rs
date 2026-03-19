pub mod aggregator;
pub mod transcript;

use crate::hash_to_curve::PairingWithG2Map;
use crate::pvss::SecretSharingWithWitness;
use crate::{pvss, VerifiedSharing};
use ark_ec::hashing::curve_maps::wb::{WBConfig, WBMap};
use ark_ec::hashing::map_to_curve_hasher::MapToCurve;
use ark_ec::pairing::Pairing;
use ark_ec::CurveGroup;
use ark_std::rand::Rng;
use ark_std::UniformRand;
use hashbrown::{HashMap, HashSet};
use transcript::{ContributionReceipt, Transcript};

/// Fully defines a protocol instance.
#[derive(Clone, Debug)]
pub struct Dkg<C: Pairing> {
    /// PVSS parameters, including the list of signers.
    pub pvss: pvss::Params<C>,
    /// Authorized dealers.
    pub dealer_pks: HashSet<C::G1Affine>,
    /// Verifies `SecretSharing`s produced by PVSS.
    pub verifier: pvss::Verifier<C>,
    /// This many of the authorized dealers has to participate to produce a "secure" transcript.
    /// With 2/3-honesty assumption, t_dkg = (1/3 + 1) of dealers
    pub t_dkg: usize,
}

/// A dealer not interested in further participation in the protocol (aggregating transcripts) can call this.
pub fn deal_and_sign<C: PairingWithG2Map, R: Rng>(
    pvss: &pvss::Params<C>,
    rng: &mut R,
    dealer: (C::ScalarField, C::G1Affine),
) -> Result<Transcript<C>, ()>
{
    let ssk = C::ScalarField::rand(rng);
    deal_and_sign_ssk(ssk, pvss, rng, dealer)
}

/// A dealer not interested in further participation in the protocol (aggregating transcripts) can call this.
pub(crate) fn deal_and_sign_ssk<C: PairingWithG2Map, R: Rng>(
    ssk: C::ScalarField,
    pvss: &pvss::Params<C>,
    rng: &mut R,
    dealer: (C::ScalarField, C::G1Affine),
) -> Result<Transcript<C>, ()>
{
    let sh = C::ScalarField::rand(rng);
    let ss = pvss.deal_secrets(ssk, sh, rng)?;
    let receipt = ContributionReceipt::<C>::sign((ssk, ss.payload.c), (sh, ss.payload.h1), dealer);
    Ok(Transcript {
        agg_ss: ss,
        receipts: vec![(receipt, 1)],
    })
}

impl<C: Pairing> Dkg<C>
where
    <C::G2 as CurveGroup>::Config: WBConfig,
    WBMap<<C::G2 as CurveGroup>::Config>: MapToCurve<C::G2>,
{
    pub fn from_pvss(
        pvss: pvss::Params<C>,
        dealer_pks: Vec<C::G1Affine>,
        t_dkg: usize,
        verifier: pvss::Verifier<C>,
    ) -> Result<Self, ()> {
        if t_dkg == 0 || t_dkg > dealer_pks.len() || verifier.config != pvss.config {
            return Err(());
        }
        let dealer_pks: HashSet<_> = dealer_pks.into_iter().collect();

        Ok(Self {
            pvss,
            dealer_pks,
            verifier,
            t_dkg,
        })
    }

    /// This method creates a new `pvss::Verifier`, so use `Self::from_pvss` if you have a compatible one.
    pub fn new(
        signer_pks: Vec<C::G2Affine>,
        t_pvss: usize,
        dealer_pks: Vec<C::G1Affine>,
        t_dkg: usize,
    ) -> Result<Self, ()> {
        let pvss = pvss::Params::<C>::new(signer_pks, t_pvss)?;
        let verifier = pvss::Verifier::new(pvss.config.clone());
        Self::from_pvss(pvss, dealer_pks, t_dkg, verifier)
    }

    pub fn deal_and_sign<R: Rng>(
        &self,
        rng: &mut R,
        dealer: (C::ScalarField, C::G1Affine),
    ) -> Result<Transcript<C>, ()> {
        deal_and_sign(&self.pvss, rng, dealer)
    }

    pub fn verify<R: Rng>(&self, transcript: &Transcript<C>, rng: &mut R) -> Result<(), ()>
    where
        <C::G2 as CurveGroup>::Config: WBConfig,
        WBMap<<C::G2 as CurveGroup>::Config>: MapToCurve<C::G2>,
    {
        // TODO: it doesn't check the dealers
        for (r, _w) in transcript.receipts.iter() {
            r.verify_all_sigs()?;
        }
        transcript.check_consistency()?;
        self.verifier
            .verify(&transcript.agg_ss, &self.pvss.signer_pks, rng)?;
        Ok(())
    }

    /// Merges `transcripts` into a new (aggregate) transcript.
    /// The result doesn't have duplicate receipts, or `0` weights.
    pub fn aggregate(transcripts: Vec<Transcript<C>>) -> Transcript<C> {
        let (ss, receipts): (Vec<_>, Vec<_>) = transcripts
            .into_iter()
            .map(|t| (t.agg_ss, t.receipts))
            .collect();
        let agg_ss = SecretSharingWithWitness::aggregate(&ss);
        let mut receipts_map: HashMap<ContributionReceipt<C>, u32> = HashMap::new();
        receipts
            .into_iter()
            .flatten()
            .filter(|(_r, w)| *w > 0)
            .for_each(|(r, w)| *receipts_map.entry(r).or_insert(0) += w);
        let receipts: Vec<_> = receipts_map.into_iter().collect();
        Transcript { agg_ss, receipts }
    }

    fn authorized_contributions(&self, t: &Transcript<C>) -> HashSet<C::G1Affine> {
        t.receipts
            .iter()
            .filter_map(|r| {
                self.dealer_pks
                    .contains(&r.0.dealer_pk)
                    .then_some(r.0.dealer_pk)
            })
            .collect()
    }

    fn enough_contributions(&self, t: &Transcript<C>) -> bool {
        self.authorized_contributions(t).len() >= self.t_dkg
    }

    pub fn finalize<R: Rng>(self, t: Transcript<C>, rng: &mut R) -> Result<VerifiedSharing<C>, ()> {
        if !self.enough_contributions(&t) {
            return Err(());
        }
        self.verify(&t, rng)?;
        Ok(VerifiedSharing {
            secret_sharing: t.agg_ss.payload,
            params: self.pvss,
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::bls::vanilla::BlsSigner;
    use crate::dkg::Dkg;
    use crate::BlsDkg;
    use ark_bls12_381::{Bls12_381, G1Affine, G2Affine};
    use ark_std::{test_rng, UniformRand};

    #[test]
    fn aggregation() {
        let rng = &mut test_rng();

        let (n, t) = (10, 7);

        let dealers: Vec<_> = (0..3).map(|_| BlsSigner::<Bls12_381>::new(rng)).collect();
        let dealer_pks: Vec<G1Affine> = dealers.iter().map(|d| d.pk_g1).collect();
        let signers_pks: Vec<_> = (0..n).map(|_| G2Affine::rand(rng)).collect();

        let dkg =
            Dkg::<Bls12_381>::new(signers_pks, t, dealer_pks.clone(), dealer_pks.len()).unwrap();
        let ss1 = dkg.deal_and_sign(rng, dealers[0].as_tuple()).unwrap();
        let ss2 = dkg.deal_and_sign(rng, dealers[1].as_tuple()).unwrap();
        let agg_ss = BlsDkg::aggregate(vec![ss1.clone(), ss1, ss2]);
        assert_eq!(agg_ss.receipts.len(), 2);
        // assert_eq!(agg_ss.receipts[0].1, 2); //TODO
    }
}
