use crate::bls::threshold::{AggThresholdSig, ThresholdVk};
use crate::bls::vanilla::BlsSigner;
use crate::bls::{dec, enc};
use ark_ec::hashing::curve_maps::wb::{WBConfig, WBMap};
use ark_ec::hashing::map_to_curve_hasher::MapToCurve;
use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ec::CurveGroup;
use ark_ff::Zero;
use ark_std::rand::Rng;
use ark_std::UniformRand;

/// Ephemeral public encryption key.
pub struct EncPk<C: Pairing> {
    pub epk_bgpk: C::G1Affine,
    pub epk_sig: C::G2Affine,
    pub epk_apk: C::G1Affine,
}

impl<C: Pairing> ThresholdVk<C>
where
    <C::G1 as CurveGroup>::Config: WBConfig,
    WBMap<<C::G1 as CurveGroup>::Config>: MapToCurve<C::G1>,
{
    pub fn initiate_key_exchange<R: Rng>(&self, id: &[u8], rng: &mut R) -> (PairingOutput<C>, EncPk<C>) {
        // esk -- ephemeral secret key
        let (esk_apk, esk_bls) = (C::ScalarField::rand(rng), C::ScalarField::rand(rng));
        // shared secret
        let ss = C::pairing(self.c * esk_apk, self.g2);
        // epk -- ephemeral public key
        let epk_bgpk = self.g1 * esk_apk;
        let epk_sig = self.g2 * esk_bls;
        let esk_bls_x_id = BlsSigner::<C>::with_sk(esk_bls).hash_and_sign(id); // TODO: computes pks
        let epk_apk = -(self.h1 * esk_apk + esk_bls_x_id.sig);
        let epk = EncPk {
            epk_bgpk: epk_bgpk.into_affine(),
            epk_sig: epk_sig.into_affine(),
            epk_apk: epk_apk.into_affine(),
        };
        (ss, epk)
    }

    pub fn check_epk(&self, epk: &EncPk<C>, id: &[u8]) -> bool {
        let id = BlsSigner::<C>::hash_to_g1(id);
        C::multi_pairing(
            &[epk.epk_bgpk, id, epk.epk_apk],
            &[self.h2, epk.epk_sig, self.g2.into_affine()],
        ).is_zero()
    }

    // TODO: unlucky method signature
    pub fn encrypt<R: Rng>(&self, id: &[u8], pt: &[u8], rng: &mut R) -> (Vec<u8>, EncPk<C>) {
        let (ss, epk) = self.initiate_key_exchange(id, rng);
        let cc = enc(pt, &ss);
        (cc, epk)
    }

    pub fn decrypt(&self, cc: &[u8], epk: &EncPk<C>, sig: &AggThresholdSig<C>) -> Vec<u8> {
        let ss = epk.complete_key_exchange(sig);
        dec(cc, &ss)
    }
}

impl<C: Pairing> EncPk<C> {
    pub fn complete_key_exchange(&self, sig: &AggThresholdSig<C>) -> PairingOutput<C> {
        // TODO: check the sig?
        C::multi_pairing(
            &[self.epk_bgpk, sig.bls_sig_with_pk.sig, self.epk_apk],
            &[sig.bgpk, self.epk_sig, sig.bls_sig_with_pk.pk],
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::bls::vanilla::BlsSigner;
    use ark_bls12_381::Bls12_381;
    use ark_ec::{AffineRepr, CurveGroup};
    use ark_std::test_rng;

    #[test]
    fn test_ibe_key_exchange() {
        let rng = &mut test_rng();

        let (n, t) = (7, 5);
        let signers: Vec<BlsSigner<Bls12_381>> = (0..n).map(|_| BlsSigner::new(rng)).collect();
        let signers_pks: Vec<_> = signers.iter().map(|s| s.bls_pk_g2).collect();
        let (tvk, sig_aggregator) = crate::tests::simulate_pvss(signers_pks, t, rng);

        let (ss, epk) = tvk.initiate_key_exchange(b"id", rng);
        assert!(tvk.check_epk(&epk, b"id"));
        assert!(!tvk.check_epk(&epk, b"id2"));

        let id_hash = BlsSigner::<Bls12_381>::hash_to_g1("id".as_bytes()).into_group();
        let sigs: Vec<_> = signers.iter().map(|s| s.sign_g1(id_hash)).collect();
        let agg_sig =  sig_aggregator.check_then_aggregate(id_hash.into_affine(), sigs.clone());
        let ss_ = epk.complete_key_exchange(&agg_sig);
        assert_eq!(ss, ss_);

        let id_hash = BlsSigner::<Bls12_381>::hash_to_g1("id2".as_bytes()).into_group();
        let sigs: Vec<_> = signers.iter().map(|s| s.sign_g1(id_hash)).collect();
        let agg_sig =  sig_aggregator.check_then_aggregate(id_hash.into_affine(), sigs.clone());
        let ss_ = epk.complete_key_exchange(&agg_sig);
        assert_ne!(ss, ss_);
    }
}