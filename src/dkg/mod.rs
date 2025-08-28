pub mod aggregator;

use crate::bls::vanilla::{hash_to_curve, sign_point, verify_on_point};
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

/// Full transcript of a (running) DKG protocol.
/// Contains a secret sharing aggregated from a number of dealers,
/// and the corresponding signatures from the dealers
/// with weights/counts (number of times the dealing has been aggregated).
#[derive(Clone, Debug)] //TODO: make a map, and eq
pub struct Transcript<C: Pairing> {
    pub agg_ss: SecretSharingWithWitness<C>,
    receipts: Vec<(ContributionReceipt<C>, u32)>,
}

/// BLS proof of possession of the secrets `(ssk, sh, sk)`,
/// corresponding to the public keys `(c = ssk.g1 = f(0).g1, h1 = sh.g1, pk = sk.g1)`.
/// All `3` signatures sign the concatenation `c || h1 || pk` of the public keys.
/// The keys `ssk` and `sh` are ephemeral (used by an honest dealer once),
/// `(sk, pk)` is the dealer's long-term keypair.
/// Together the signatures show that who knows `sk`, knows also `ssk` and `sh`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ContributionReceipt<C: Pairing> {
    // BLS public keys in G1
    c: C::G1Affine,
    h1: C::G1Affine,
    dealer_pk: C::G1Affine,
    // BLS signatures in G2
    sig_c: C::G2Affine,
    sig_h1: C::G2Affine,
    sig_pk: C::G2Affine,
}

impl<C: Pairing> ContributionReceipt<C>
where
    <C::G2 as CurveGroup>::Config: WBConfig,
    WBMap<<C::G2 as CurveGroup>::Config>: MapToCurve<C::G2>,
{
    fn sign(c: (C::ScalarField, C::G1Affine), h1: (C::ScalarField, C::G1Affine), dealer: (C::ScalarField, C::G1Affine)) -> Self {
        let public_keys = (c.1, h1.1, dealer.1);
        let message_hash = hash_to_curve::<C::G2, _>(public_keys);
        Self {
            c: c.1,
            h1: h1.1,
            dealer_pk: dealer.1,
            sig_c: sign_point(c.0, message_hash),
            sig_h1: sign_point(h1.0, message_hash),
            sig_pk: sign_point(dealer.0, message_hash),
        }
    }

    fn hash_pks(&self) -> C::G2Affine {
        let public_keys = (self.c, self.h1, self.dealer_pk);
        hash_to_curve::<C::G2, _>(public_keys)
    }

    fn verify_all_sigs(&self) -> Result<(), ()> {
        let message_hash = self.hash_pks();
        let g1 = C::G1::generator();
        if verify_on_point::<C>(self.sig_c, message_hash, self.c, g1)
            && verify_on_point::<C>(self.sig_h1, message_hash, self.h1, g1)
            && verify_on_point::<C>(self.sig_pk, message_hash, self.dealer_pk, g1) {
            Ok(())
        } else {
            Err(())
        }
    }

    // fn verify_dealer(&self) -> bool {
    //     let message_hash = self.hash_pks();
    //     let g1 = C::G1::generator();
    //     verify_on_point::<C>(self.sig_pk, message_hash, self.dealer_pk, g1)
    // }
}

impl<C: Pairing> Hash for ContributionReceipt<C> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.c.hash(state);
        self.h1.hash(state);
        self.dealer_pk.hash(state);
        self.sig_c.hash(state);
        self.sig_h1.hash(state);
        self.sig_pk.hash(state);
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Dkg<C: Pairing> {
    pub pvss: pvss::Params<C>,
    pub dealer_pks: HashSet<C::G1Affine>,
    pub t_dkg: usize,
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
        let ssk = C::ScalarField::rand(rng);
        let sh = C::ScalarField::rand(rng);
        let pvss = self.pvss.deal_secrets(ssk, sh, rng);
        let receipt = ContributionReceipt::sign(
            (ssk, pvss.payload.c),
            (sh, pvss.payload.h1),
            dealer,
        );
        Transcript {
            agg_ss: pvss,
            receipts: vec![(receipt, 1)],
        }
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

    pub fn finalize<R: Rng>(self, t: Transcript<C>, rng: &mut R) -> Result<ThresholdCrypto<C>,()> {
        let v = pvss::Verifier::new(self.pvss.config.clone());
        self.verify(&t, &v, rng)?;
        if !self.enough_dealers(&t) {
            return Err(())
        }
        Ok(ThresholdCrypto {
            secret_sharing: t.agg_ss.payload,
            params: self.pvss,
        })
    }
}

impl<C: Pairing> Transcript<C> {
    pub fn list_dealers(&self) -> Vec<C::G1Affine> {
        self.receipts.iter()
            .map(|(r, _w)| r.dealer_pk)
            .collect()
    }

    /// `c`s and `h1`s in the receipts sum up to `c` and `h1` in the secret sharing.
    pub fn check_consistency(&self) -> Result<(), ()> {
        let cs: Vec<_> = self.receipts.iter().map(|(r, _w)| r.c).collect();
        let h1s: Vec<_> = self.receipts.iter().map(|(r, _w)| r.h1).collect();
        let ws: Vec<_> = self.receipts.iter().map(|(_, w)| C::ScalarField::from(*w)).collect();
        let c = C::G1::msm(&cs, &ws).unwrap();
        if c.into_affine() != self.agg_ss.payload.c {
            return Err(());
        }
        let h1 = C::G1::msm(&h1s, &ws).unwrap();
        if h1.into_affine() != self.agg_ss.payload.h1 {
            return Err(());
        }
        Ok(())
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