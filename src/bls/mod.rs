use ark_ec::pairing::{Pairing, PairingOutput};
use ark_serialize::CanonicalSerialize;
use sha2::{Digest, Sha256};

/// Threshold BLS verification and the VUF.
pub mod threshold;

/// Vanilla BLS signature scheme.
pub mod vanilla;
pub mod ibe;
pub mod pke;


pub fn enc<C: Pairing>(pt: &[u8], ss: &PairingOutput<C>) -> Vec<u8> {
    let mut key = vec![0; ss.compressed_size()];
    ss.serialize_compressed(&mut key[..]).unwrap();
    let mut hasher = Sha256::new();
    hasher.update(key);
    let mut hash = hasher.finalize().to_vec();
    hash.iter_mut()
        .zip(pt.iter())
        .for_each(|(x1, &x2)| *x1 ^= x2);
    hash
}

pub fn dec<C: Pairing>(cc: &[u8], ss: &PairingOutput<C>) -> Vec<u8> {
    enc(cc, ss)
}

#[cfg(test)]
mod tests {
    use crate::bls::{dec, enc};
    use ark_bls12_381::Bls12_381;
    use ark_ec::pairing::PairingOutput;
    use ark_std::rand::Rng;
    use ark_std::{test_rng, UniformRand};

    #[test]
    fn test_xor() {
        let rng = &mut test_rng();

        let ss = PairingOutput::<Bls12_381>::rand(rng);

        let pt = pt(rng);

        let cc = enc(&pt, &ss);
        let pt_ = dec(&cc, &ss);

        assert_eq!(pt_, pt);
    }

    pub fn pt<R: Rng>(rng: &mut R) -> Vec<u8> {
        let mut pt = vec![0u8; 32];
        rng.fill(&mut pt[..]);
        pt
    }
}