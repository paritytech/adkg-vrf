use crate::bls::enc;
use crate::bls::ibe::EncPk;
use crate::bls::threshold::ThresholdVk;
use crate::bls::vanilla::{BlsSigner, StandaloneSig};
use ark_ec::hashing::curve_maps::wb::{WBConfig, WBMap};
use ark_ec::hashing::map_to_curve_hasher::MapToCurve;
use ark_ec::pairing::Pairing;
use ark_ec::{AffineRepr, CurveGroup};
use ark_std::rand::Rng;
use ark_std::UniformRand;

impl<C: Pairing> ThresholdVk<C>
where
    <C::G1 as CurveGroup>::Config: WBConfig,
    WBMap<<C::G1 as CurveGroup>::Config>: MapToCurve<C::G1>,
{
    // Doesn't call Self::initiate_key_exchange as here the identity depends on the ciphertext and parts of the epk.
    pub fn encrypt_to_threshold<R: Rng>(&self, pt: &[u8], rng: &mut R) -> (Vec<u8>, EncPk<C>) {
        // esk -- ephemeral secret key
        let (esk_apk, esk_bls) = (C::ScalarField::rand(rng), C::ScalarField::rand(rng));
        // shared secret
        let ss = C::pairing(self.c * esk_apk, self.g2);
        // epk -- ephemeral public key
        let epk_bgpk = self.g1 * esk_apk;
        let epk_sig = self.g2 * esk_bls;
        // ciphertext
        let cc = enc(pt, &ss);
        // identity to encrypt to
        let id = (cc.clone(), epk_bgpk, epk_sig);
        let esk_bls_x_id = BlsSigner::<C>::with_sk(esk_bls).hash_and_sign(id); // TODO: computes pks
        let epk_apk = -(self.h1 * esk_apk + esk_bls_x_id.sig);
        let epk = EncPk {
            epk_bgpk: epk_bgpk.into_affine(),
            epk_sig: epk_sig.into_affine(),
            epk_apk: epk_apk.into_affine(),
        };
        (cc, epk)
    }

    // Signs the id
    pub fn partial_decrypt(&self, sk: C::ScalarField, cc: &[u8], epk: &EncPk<C>) -> StandaloneSig<C> {
        let id = (cc.clone(), epk.epk_bgpk, epk.epk_sig);
        let id_hash = BlsSigner::<C>::hash_to_g1(id);
        BlsSigner::<C>::with_sk(sk).sign_g1(id_hash.into_group())
    }
}

#[cfg(test)]
mod tests {
    use crate::bls::threshold::AggThresholdSig;
    use crate::bls::vanilla::BlsSigner;
    use ark_bls12_381::Bls12_381;
    use ark_std::test_rng;

    #[test]
    fn test_threshold_pke() {
        let rng = &mut test_rng();

        let (n, t) = (1, 1);
        let signers: Vec<BlsSigner<Bls12_381>> = (0..n).map(|_| BlsSigner::new(rng)).collect();
        let signers_pks: Vec<_> = signers.iter().map(|s| s.bls_pk_g2).collect();
        let (tvk, bgpk) = crate::tests::simulate_pvss::<Bls12_381, _>(signers_pks, t, rng);

        let signer = &signers[0];
        let bgpk = bgpk[0];

        let pt = crate::bls::tests::pt(rng);

        let (cc, epk) = tvk.encrypt_to_threshold(&pt, rng);

        let sig = tvk.partial_decrypt(signer.sk, &cc, &epk);
        let sig = AggThresholdSig { bls_sig_with_pk: sig, bgpk };

        let pt_ = tvk.decrypt(&cc, &epk, &sig);
        assert_eq!(pt, pt_);
    }
}