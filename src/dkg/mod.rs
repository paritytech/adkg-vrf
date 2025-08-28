pub mod aggregator;
pub mod transcript;

use crate::pvss::SecretSharingWithWitness;
use crate::{pvss, ThresholdCrypto};
use ark_ec::hashing::curve_maps::wb::{WBConfig, WBMap};
use ark_ec::hashing::map_to_curve_hasher::MapToCurve;
use ark_ec::pairing::Pairing;
use ark_ec::{CurveGroup, PrimeGroup, VariableBaseMSM};
use ark_std::rand::Rng;
use ark_std::UniformRand;
use hashbrown::{HashMap, HashSet};
use std::hash::{Hash, Hasher};
use transcript::{ContributionReceipt, Transcript};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Dkg<C: Pairing> {
    pub pvss: pvss::Params<C>,
    pub dealer_pks: HashSet<C::G1Affine>,
    pub t_dkg: usize,
}

/// A dealer not interested in further participation in the protocol (aggregating transcripts) can call this.
pub fn deal_and_sign<C: Pairing, R: Rng>(pvss: &pvss::Params<C>, rng: &mut R, dealer: (C::ScalarField, C::G1Affine)) -> Transcript<C>
where
    <C::G2 as CurveGroup>::Config: WBConfig,
    WBMap<<C::G2 as CurveGroup>::Config>: MapToCurve<C::G2>,
{
    let ssk = C::ScalarField::rand(rng);
    let sh = C::ScalarField::rand(rng);
    let ss = pvss.deal_secrets(ssk, sh, rng);
    let receipt = ContributionReceipt::<C>::sign(
        (ssk, ss.payload.c),
        (sh, ss.payload.h1),
        dealer,
    );
    Transcript {
        agg_ss: ss,
        receipts: vec![(receipt, 1)],
    }
}

impl<C: Pairing> Dkg<C>
where
    <C::G2 as CurveGroup>::Config: WBConfig,
    WBMap<<C::G2 as CurveGroup>::Config>: MapToCurve<C::G2>,
{
    pub fn from_pvss(pvss: pvss::Params<C>, dealer_pks: Vec<C::G1Affine>, t_dkg: usize) -> Result<Self, ()> {
        let dealer_pks: HashSet<_> = dealer_pks.into_iter().collect();
        if t_dkg == 0 || t_dkg > dealer_pks.len() {
            return Err(());
        }
        Ok(Self { pvss, dealer_pks, t_dkg })
    }

    pub fn new(signer_pks: Vec<C::G2Affine>, t_pvss: usize, dealer_pks: Vec<C::G1Affine>, t_dkg: usize) -> Result<Self, ()> {
        let pvss = pvss::Params::<C>::new(signer_pks, t_pvss)?;
        Self::from_pvss(pvss, dealer_pks, t_dkg)
    }

    pub fn deal_and_sign<R: Rng>(&self, rng: &mut R, dealer: (C::ScalarField, C::G1Affine)) -> Transcript<C> {
        deal_and_sign(&self.pvss, rng, dealer)
    }

    pub fn verify<R: Rng>(
        &self,
        transcript: &Transcript<C>,
        pvss_verifier: &pvss::Verifier<C>,
        rng: &mut R,
    ) -> Result<(), ()>
    where
        <C::G2 as CurveGroup>::Config: WBConfig,
        WBMap<<C::G2 as CurveGroup>::Config>: MapToCurve<C::G2>,
    {
        // TODO: it doesn't check the dealers
        for (r, _w) in transcript.receipts.iter() {
            r.verify_all_sigs()?;
        }
        transcript.check_consistency()?;
        pvss_verifier.verify(&transcript.agg_ss, &self.pvss.signer_pks, rng)?;
        Ok(())
    }

    pub fn aggregate(transcripts: Vec<Transcript<C>>) -> Transcript<C> {
        let (pvss, witness): (Vec<_>, Vec<_>) = transcripts.into_iter()
            .map(|t| (t.agg_ss, t.receipts))
            .collect();
        let agg_pvss = SecretSharingWithWitness::aggregate(&pvss);
        let mut weights: HashMap<ContributionReceipt<C>, u32> = HashMap::new();
        witness.into_iter()
            .flatten()
            .for_each(|(c, w)| *weights.entry(c).or_insert(0) += w);
        let receipts: Vec<_> = weights.into_iter().collect();
        Transcript { agg_ss: agg_pvss, receipts }
    }

    fn contributed_dealers(&self, t: &Transcript<C>) -> HashSet<C::G1Affine> {
        t.receipts.iter()
            .filter_map(|r| self.dealer_pks.contains(&r.0.dealer_pk).then_some(r.0.dealer_pk))
            .collect()
    }

    fn enough_dealers(&self, t: &Transcript<C>) -> bool {
        self.contributed_dealers(t).len() >= self.t_dkg
    }

    pub fn finalize<R: Rng>(self, t: Transcript<C>, rng: &mut R) -> Result<ThresholdCrypto<C>, ()> {
        let v = pvss::Verifier::new(self.pvss.config.clone());
        self.verify(&t, &v, rng)?;
        if !self.enough_dealers(&t) {
            return Err(());
        }
        Ok(ThresholdCrypto {
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

        let dealers: Vec<_> = (0..3)
            .map(|_| BlsSigner::<Bls12_381>::new(rng))
            .collect();
        let dealer_pks: Vec<G1Affine> = dealers.iter()
            .map(|d| d.bls_pk_g1)
            .collect();
        let signers_pks: Vec<_> = (0..n)
            .map(|_| G2Affine::rand(rng))
            .collect();

        let dkg = Dkg::<Bls12_381>::new(signers_pks, t, dealer_pks.clone(), dealer_pks.len()).unwrap();
        let ss1 = dkg.deal_and_sign(rng, dealers[0].pk_in_g1());
        let ss2 = dkg.deal_and_sign(rng, dealers[1].pk_in_g1());
        let agg_ss = BlsDkg::aggregate(vec![ss1.clone(), ss1, ss2]);
        assert_eq!(agg_ss.receipts.len(), 2);
        // assert_eq!(agg_ss.receipts[0].1, 2); //TODO
    }
}