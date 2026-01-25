use ark_ec::CurveGroup;
use ark_ff::One;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::Rng;
use ark_std::UniformRand;
use ark_std::{iter, vec, vec::Vec};

use crate::utils::powers;

/// Batch-verifiable proofs of knowledge of discrete logarithms
/// of a number of `points` with respect to the same `base`.
#[derive(CanonicalSerialize)]
pub struct Instance<G: CurveGroup> {
    pub base: G,
    pub points: Vec<G>,
}

pub struct Statement<G: CurveGroup> {
    pub instance: Instance<G>,
    pub dlogs: Vec<G::ScalarField>,
}

#[derive(Clone, Debug, PartialEq, Eq, CanonicalDeserialize, CanonicalSerialize)]
pub struct Proof<G: CurveGroup> {
    r: G,
    s: G::ScalarField,
}

fn c<G: CurveGroup>(instance: &Instance<G>, r: &G) -> G::ScalarField {
    let mut t = ark_transcript::Transcript::new_blank();
    t.append(instance);
    t.append(r);
    t.challenge(b"whatever").read_reduce()
}

impl<G: CurveGroup> Statement<G> {
    pub fn prove<R: Rng>(&self, rng: &mut R) -> Proof<G> {
        let r = G::ScalarField::rand(rng);
        let r_big = self.instance.base * r;
        let c = c(&self.instance, &r_big);
        let s: G::ScalarField = self
            .dlogs
            .iter()
            .zip(powers(c).skip(1))
            .map(|(exp, c)| c * exp)
            .sum();
        let s = r - s;
        Proof { r: r_big, s }
    }
}

impl<G: CurveGroup> Instance<G> {
    #[must_use]
    pub fn verify(&self, proof: &Proof<G>) -> bool {
        let c = c(self, &proof.r);
        let p: G = self
            .points
            .iter()
            .zip(powers(c).skip(1))
            .map(|(&r, ci)| r * ci)
            .sum();
        proof.r == self.base * proof.s + p
    }

    // TODO: check https://eprint.iacr.org/2022/222.pdf
    #[must_use]
    pub fn batch_verify<R: Rng>(claims: &[(Instance<G>, Proof<G>)], rng: &mut R) -> bool {
        let l = G::ScalarField::rand(rng);
        let coeffs: Vec<_> = claims
            .iter()
            .flat_map(|(x, pi)| {
                let mut tuple = vec![pi.s, -G::ScalarField::one()];
                tuple.extend(powers(c(x, &pi.r)).skip(1).take(x.points.len()));
                tuple
            })
            .collect();
        let ls: Vec<_> = claims
            .iter()
            .zip(powers(l))
            .flat_map(|((x, _), li)| iter::repeat(li).take(x.points.len() + 2))
            .collect();
        assert_eq!(coeffs.len(), ls.len());
        let coeffs: Vec<_> = coeffs.into_iter().zip(ls).map(|(ci, li)| ci * li).collect();
        let bases: Vec<_> = claims
            .iter()
            .flat_map(|(x, pi)| [&vec![x.base, pi.r], x.points.as_slice()].concat())
            .collect();
        let bases = G::normalize_batch(&bases);
        G::msm(&bases, &coeffs).unwrap().is_zero()
    }
}

#[cfg(test)]
mod tests {
    use ark_ec::PrimeGroup;
    use ark_std::vec::Vec;
    use ark_std::{test_rng, UniformRand};

    use crate::koe::{Instance, Statement};

    #[test]
    fn koe() {
        let rng = &mut test_rng();

        let k = 10;
        let base = ark_bls12_381::G1Projective::generator();
        let exponents = (0..k)
            .map(|_| ark_bls12_381::Fr::rand(rng))
            .collect::<Vec<_>>();
        let results = exponents.iter().map(|exp| base * exp).collect::<Vec<_>>();

        let instance = Instance {
            base,
            points: results,
        };
        let statement = Statement {
            instance,
            dlogs: exponents,
        };

        let proof = statement.prove(rng);
        assert!(statement.instance.verify(&proof));

        assert!(Instance::batch_verify(&[(statement.instance, proof)], rng));
    }
}
