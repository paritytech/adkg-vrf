use ark_ec::hashing::curve_maps::wb::{WBConfig, WBMap};
use ark_ec::hashing::map_to_curve_hasher::{MapToCurve, MapToCurveBasedHasher};
use ark_ec::hashing::{HashToCurve, HashToCurveError};
use ark_ec::pairing::Pairing;
use ark_ec::CurveGroup;
use ark_ff::field_hashers::DefaultFieldHasher;
use ark_serialize::CanonicalSerialize;
use sha2::Sha256;

pub trait CurveWithPairingAndHash: PairingWithG1Map + PairingWithG2Map {}

pub trait PairingWithG1Map: Pairing {
    const DOMAIN: &'static [u8];
    fn hash_to_g1(bytes: &[u8]) -> Result<Self::G1Affine, HashToCurveError>;

    fn hash_serializable_to_g1<M: CanonicalSerialize>(msg: &M) -> Result<Self::G1Affine, HashToCurveError> {
        let mut bytes = vec![0; msg.compressed_size()];
        msg.serialize_compressed(&mut bytes[..]).unwrap();
        Self::hash_to_g1(&bytes)
    }
}

pub trait PairingWithG2Map: Pairing {
    const DOMAIN: &'static [u8];
    fn hash_to_g2(bytes: &[u8]) -> Result<Self::G2Affine, HashToCurveError>;

    fn hash_serializable_to_g2<M: CanonicalSerialize>(msg: &M) -> Result<Self::G2Affine, HashToCurveError> {
        let mut bytes = vec![0; msg.compressed_size()];
        msg.serialize_compressed(&mut bytes[..]).unwrap();
        Self::hash_to_g2(&bytes)
    }
}

type FieldHasher = DefaultFieldHasher<Sha256, 128>;

type G1Config<C> = <<C as Pairing>::G1 as CurveGroup>::Config;
type WbToG1<C> = WBMap<G1Config<C>>;

impl<C> PairingWithG1Map for C
where
    C: Pairing,
    G1Config<C>: WBConfig,
    WbToG1<C>: MapToCurve<C::G1>,
{
    const DOMAIN: &'static [u8] = b"adkg-vrf-hash-t-g1";

    // TODO: there shouldn't be any errors afaik
    fn hash_to_g1(bytes: &[u8]) -> Result<C::G1Affine, HashToCurveError> {
        let wb_hasher = MapToCurveBasedHasher::<C::G1, FieldHasher, WbToG1<C>>::new(Self::DOMAIN)?;
        wb_hasher.hash(&bytes)
    }
}

type G2Config<C> = <<C as Pairing>::G2 as CurveGroup>::Config;
type WbToG2<C> = WBMap<G2Config<C>>;

impl<C: Pairing> PairingWithG2Map for C
where
    G2Config<C>: WBConfig,
    WbToG2<C>: MapToCurve<C::G2>,
{
    const DOMAIN: &'static [u8] = b"adkg-vrf-hash-t-g2";

    // TODO: there shouldn't be any errors afaik
    fn hash_to_g2(bytes: &[u8]) -> Result<C::G2Affine, HashToCurveError> {
        let wb_hasher = MapToCurveBasedHasher::<C::G2, FieldHasher, WbToG2<C>>::new(Self::DOMAIN)?;
        wb_hasher.hash(&bytes)
    }
}

impl<C> CurveWithPairingAndHash for C
where
    C: PairingWithG1Map,
    C: PairingWithG2Map,
{}

#[cfg(test)]
mod tests {
    use crate::hash_to_curve::CurveWithPairingAndHash;

    fn test_g1_map<C: CurveWithPairingAndHash>() {
        assert!(C::hash_to_g1(b"map-to-g1").is_ok());
        assert!(C::hash_to_g2(b"map-to-g2").is_ok());
    }

    #[test]
    fn test_bls12_381_g1_map() {
        test_g1_map::<ark_bls12_381::Bls12_381>();
    }
}