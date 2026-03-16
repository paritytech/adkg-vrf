use ark_ec::pairing::Pairing;
use ark_ec::{AffineRepr, CurveGroup, PrimeGroup};
use ark_ff::Zero;
use ark_serialize::CanonicalSerialize;
use ark_std::rand::Rng;
use ark_std::UniformRand;

use crate::PairingWithG1Map;

pub fn sign_point<C: AffineRepr>(sk: C::ScalarField, point: C) -> C {
    let sig = point * sk;
    sig.into_affine()
}

pub fn verify_on_point<C: Pairing>(
    sig: C::G2Affine,
    point: C::G2Affine,
    pk: C::G1Affine,
    g1: C::G1,
) -> bool {
    let minus_g1 = (-g1).into_affine();
    C::multi_pairing([minus_g1, pk], [sig, point]).is_zero()
}

#[derive(Clone)]
pub struct BlsSigner<C: Pairing> {
    pub sk: C::ScalarField,
    pub bls_pk_g1: C::G1Affine,
    pub bls_pk_g2: C::G2Affine,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StandaloneSig<C: Pairing> {
    pub sig: C::G1Affine,
    pub pk: C::G2Affine,
}

impl<C: Pairing> BlsSigner<C> {
    pub fn new<R: Rng>(rng: &mut R) -> Self {
        let sk = C::ScalarField::rand(rng);
        Self::with_sk(sk)
    }

    pub fn with_sk(sk: C::ScalarField) -> Self {
        let bls_pk_g1 = C::G1::generator() * sk;
        let bls_pk_g2 = C::G2::generator() * sk;
        let bls_pk_g1 = bls_pk_g1.into_affine();
        let bls_pk_g2 = bls_pk_g2.into_affine();
        Self {
            sk,
            bls_pk_g1,
            bls_pk_g2,
        }
    }

    pub fn sign_g1(&self, m: C::G1) -> StandaloneSig<C> {
        let sig = m * self.sk;
        let sig = sig.into_affine();
        StandaloneSig {
            sig,
            pk: self.bls_pk_g2,
        }
    }

    pub fn sign_g2(&self, m: C::G2) -> C::G2Affine {
        let sig = m * self.sk;
        let sig = sig.into_affine();
        sig
    }

    pub fn pk_in_g1(&self) -> (C::ScalarField, C::G1Affine) {
        (self.sk, self.bls_pk_g1)
    }
}

impl<C: PairingWithG1Map> BlsSigner<C>
{
    pub fn hash_to_g1<M: CanonicalSerialize>(msg: M) -> C::G1Affine {
        C::hash_serializable_to_g1(&msg).unwrap()
    }

    pub fn hash_and_sign<M: CanonicalSerialize>(&self, m: M) -> StandaloneSig<C> {
        let m_hash_g1 = Self::hash_to_g1(m);
        self.sign_g1(m_hash_g1.into_group())
    }
}

impl<C: Pairing> StandaloneSig<C> {
    pub fn verify_unoptimized(&self, m: C::G1, g2: C::G2Affine) {
        assert_eq!(C::pairing(self.sig, g2), C::pairing(m.into(), self.pk));
    }
}
