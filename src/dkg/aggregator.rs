use crate::dkg;
use crate::dkg::{ContributionReceipt, Transcript};
use crate::pvss::{SecretSharingWithWitness, Verifier};
use ark_ec::hashing::curve_maps::wb::{WBConfig, WBMap};
use ark_ec::hashing::map_to_curve_hasher::MapToCurve;
use ark_ec::pairing::Pairing;
use ark_ec::CurveGroup;
use ark_std::rand::Rng;
use hashbrown::{HashMap, HashSet};

/// Aggregates DKG transcripts into a fully-verifiable DKG transcript.
#[derive(Clone)]
pub struct TranscriptAggregator<C: Pairing> {
    dkg: dkg::Dkg<C>,
    dealer_pks: HashSet<C::G1Affine>,

    agg_ss: Option<SecretSharingWithWitness<C>>,
    receipts: HashMap<ContributionReceipt<C>, u32>,
    pvss_verifier: Verifier<C>,
}

impl<C: Pairing> TranscriptAggregator<C>
where
    <C::G2 as CurveGroup>::Config: WBConfig,
    WBMap<<C::G2 as CurveGroup>::Config>: MapToCurve<C::G2>,
{
    pub fn new(dkg: dkg::Dkg<C>, dealer_pks: Vec<C::G1Affine>) -> Self {
        let pvss_verifier = Verifier::new(dkg.pvss.config.clone());
        Self {
            dkg,
            dealer_pks: dealer_pks.iter().copied().collect(),
            agg_ss: None,
            receipts: Default::default(),
            pvss_verifier,
        }
    }

    pub fn add<R: Rng>(&mut self, transcript: Transcript<C>, rng: &mut R) -> Result<(), ()> {
        let new_receipts: Vec<_> = transcript.receipts.iter()
            .filter(|(_r, w)| *w > 0)
            .collect();
        let aggregated_pks: HashSet<C::G1Affine> = self.receipts.keys()
            .map(|r| r.dealer_pk)
            .collect();
        let missing_pks: HashSet<C::G1Affine> = self.dealer_pks.difference(&aggregated_pks)
            .copied()
            .collect();
        let has_new_pks = new_receipts.iter()
            .any(|(r, _w)| missing_pks.contains(&r.dealer_pk));
        if !has_new_pks {
            return Ok(());
        }

        let receipts_to_verify: HashSet<ContributionReceipt<C>> = new_receipts.iter()
            .filter_map(|(r, _w)| (!self.receipts.contains_key(r)).then_some(r.clone()))
            .collect();

        for r in receipts_to_verify.iter() {
            r.verify_all_sigs()?;
        }

        transcript.check_consistency()?;
        self.pvss_verifier.verify(&transcript.agg_ss, &self.dkg.pvss.signer_pks, rng)?;

        for (r, w) in new_receipts {
            *self.receipts.entry(r.clone()).or_insert(0) += w;
        }

        if self.agg_ss.is_none() {
            self.agg_ss = Some(transcript.agg_ss);
        } else {
            let new_agg_ss = self.agg_ss.clone().unwrap().aggregate_with(vec![transcript.agg_ss]);
            self.agg_ss = Some(new_agg_ss);
        }
        Ok(())
    }

    pub fn get_transcript(&self) -> Transcript<C> {
        Transcript {
            agg_ss: self.agg_ss.clone().unwrap(),
            receipts: self.receipts.clone().into_iter().collect(),
        }
    }

    fn aggregated_dealer_pks(&self) -> HashSet<C::G1Affine> {
        self.receipts.keys()
            .filter_map(|r| self.dealer_pks.contains(&r.dealer_pk).then_some(r.dealer_pk))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use crate::bls::vanilla::BlsSigner;
    use crate::dkg::Dkg;
    use crate::BlsTranscriptAggregator;
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
        let ss3 = dkg.deal_and_sign(rng, dealers[2].pk_in_g1());

        let agg = BlsTranscriptAggregator::new(dkg.clone(), dealer_pks);

        // 1. aggregation of a single transcript is still the same transcript
        let mut agg1 = agg.clone();
        assert!(agg1.add(ss1.clone(), rng).is_ok());
        let expected_ss1 = agg1.get_transcript();
        assert_eq!(expected_ss1.agg_ss, ss1.agg_ss);
        assert_eq!(expected_ss1.receipts, ss1.receipts); // order of receipts may differ for n > 1

        // 2. same transcript doesn't get aggregated twice
        assert!(agg1.add(ss1.clone(), rng).is_ok());
        assert_eq!(expected_ss1.agg_ss, ss1.agg_ss);
        assert_eq!(expected_ss1.receipts, ss1.receipts); // order of receipts may differ for n > 1

        // 3. can aggregate 2 singletons
        assert!(agg1.add(ss2, rng).is_ok());
        let ss12 = agg1.get_transcript();
        assert!(dkg.verify(&ss12, &agg.pvss_verifier, rng).is_ok());
        assert_eq!(ss12.receipts.len(), 2);

        // 4. can aggregate 123 = 12 + 3
        let mut agg2 = agg.clone();
        assert!(agg2.add(ss12, rng).is_ok());
        assert!(agg2.add(ss3, rng).is_ok());
        let ss123 = agg2.get_transcript();
        assert!(dkg.verify(&ss123, &agg.pvss_verifier, rng).is_ok());
        assert_eq!(ss123.receipts.len(), 3);

        // 5. 12 + 123
        assert!(agg1.add(ss123, rng).is_ok());
        let ss = agg1.get_transcript();
        assert!(dkg.verify(&ss, &agg.pvss_verifier, rng).is_ok());
        assert_eq!(ss.receipts.len(), 3);
        assert_eq!(ss.receipts.iter().map(|(_, w)| w).sum::<u32>(), 5);  // TODO: map

    }
}