use crate::bls::vanilla::{hash_to_curve, sign_point, verify_on_point};
use crate::pvss::SecretSharingWithWitness;
use ark_ec::hashing::curve_maps::wb::{WBConfig, WBMap};
use ark_ec::hashing::map_to_curve_hasher::MapToCurve;
use ark_ec::pairing::Pairing;
use ark_ec::{CurveGroup, PrimeGroup, VariableBaseMSM};
use std::hash::{Hash, Hasher};

/// Full transcript of a (running) DKG protocol.
/// Contains a secret sharing aggregated from a number of dealers,
/// and the corresponding signatures from the dealers
/// with weights/counts (number of times the dealing has been aggregated).
#[derive(Clone, Debug)] //TODO: make a map, and eq
pub struct Transcript<C: Pairing> {
    pub agg_ss: SecretSharingWithWitness<C>,
    pub receipts: Vec<(ContributionReceipt<C>, u32)>,
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
    pub dealer_pk: C::G1Affine,
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
    pub fn sign(
        c: (C::ScalarField, C::G1Affine),
        h1: (C::ScalarField, C::G1Affine),
        dealer: (C::ScalarField, C::G1Affine),
    ) -> Self {
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

    pub fn verify_all_sigs(&self) -> Result<(), ()> {
        let message_hash = self.hash_pks();
        let g1 = C::G1::generator();
        if verify_on_point::<C>(self.sig_c, message_hash, self.c, g1)
            && verify_on_point::<C>(self.sig_h1, message_hash, self.h1, g1)
            && verify_on_point::<C>(self.sig_pk, message_hash, self.dealer_pk, g1)
        {
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

impl<C: Pairing> Transcript<C> {
    pub fn list_dealers(&self) -> Vec<C::G1Affine> {
        self.receipts.iter().map(|(r, _w)| r.dealer_pk).collect()
    }

    /// `c`s and `h1`s in the receipts sum up to `c` and `h1` in the secret sharing.
    pub fn check_consistency(&self) -> Result<(), ()> {
        let cs: Vec<_> = self.receipts.iter().map(|(r, _w)| r.c).collect();
        let h1s: Vec<_> = self.receipts.iter().map(|(r, _w)| r.h1).collect();
        let ws: Vec<_> = self
            .receipts
            .iter()
            .map(|(_, w)| C::ScalarField::from(*w))
            .collect();
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
