use crate::bls::threshold::{AggThresholdSig, ThresholdVk};
use crate::bls::vanilla::BlsSigner;
use ark_ec::hashing::curve_maps::wb::{WBConfig, WBMap};
use ark_ec::hashing::map_to_curve_hasher::MapToCurve;
use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ec::CurveGroup;
use ark_ff::Zero;
use ark_std::rand::Rng;
use ark_std::UniformRand;

/// Ephemeral public encryption key.
pub struct EncPk<C: Pairing> {
    epk_bgpk: C::G1Affine,
    epk_sig: C::G2Affine,
    epk_apk: C::G1Affine,
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
        let ss = C::pairing(self.c * esk_apk, self.g2.into());
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
    use crate::bls::threshold::{AggThresholdSig, ThresholdVk};
    use crate::bls::vanilla::BlsSigner;
    use crate::pvss;
    use ark_bls12_381::Bls12_381;
    use ark_std::test_rng;

    #[test]
    fn test_ibe_key_exchange() {
        let rng = &mut test_rng();

        let signer = BlsSigner::new(rng);

        let pvss = pvss::Params::<Bls12_381>::new(vec![signer.bls_pk_g2], 1).unwrap();
        let share = pvss.deal(rng).unwrap();
        let bgpk = share.payload.bgpk[0];
        let tvk = ThresholdVk::from_share(&share.payload);
        let (ss, epk) = tvk.initiate_key_exchange(b"id", rng);
        assert!(tvk.check_epk(&epk, b"id"));
        assert!(!tvk.check_epk(&epk, b"id2"));
        let sig = signer.hash_and_sign("id".as_bytes());
        let sig = AggThresholdSig { bls_sig_with_pk: sig, bgpk };
        let ss_ = epk.complete_key_exchange(&sig);
        assert_eq!(ss, ss_);
        let sig = signer.hash_and_sign("id2".as_bytes());
        let sig = AggThresholdSig { bls_sig_with_pk: sig, bgpk };
        let ss_ = epk.complete_key_exchange(&sig);
        assert_ne!(ss, ss_);
    }
}