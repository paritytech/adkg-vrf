use crate::pvss::{Config, SecretSharingWithWitness};
use crate::utils::BarycentricDomain;
use ark_ec::pairing::Pairing;
use ark_ec::VariableBaseMSM;
use ark_ff::{Field, One, Zero};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::rand::Rng;
use ark_std::vec::Vec;
use ark_std::{end_timer, start_timer, UniformRand};
use derivative::Derivative;

/// Precomputed barycentric weights to facilitate interpolation.
/// Depend only on `(t,n)` so can be reused between the ceremonies.
#[derive(Derivative, CanonicalSerialize, CanonicalDeserialize)]
#[derivative(Clone, Debug)]
pub struct Verifier<C: Pairing> {
    config: Config<C>,
    #[derivative(Debug="ignore")]
    domain_size_n: BarycentricDomain<C::ScalarField>,
    #[derivative(Debug="ignore")]
    domain_size_t: BarycentricDomain<C::ScalarField>,
}

impl<C: Pairing> Verifier<C> {
    /// TODO: 1. can be computed faster
    /// TODO: 2. can keep lis_at_0
    /// TODO: 3. lis_at_0 can be computed faster
    pub fn new(config: Config<C>) -> Self {
        let _t = start_timer!(|| "Interpolation");
        let domain_size_n = BarycentricDomain::of_size(config.domain, config.n);
        let domain_size_t = BarycentricDomain::of_size(config.domain, config.t);
        end_timer!(_t);
        Self {
            config,
            domain_size_n,
            domain_size_t,
        }
    }

    /// Uses Fiat-Shamir randomness to verify the DKG `transcript`, and returns the payload, if valid.
    /// Is useful for on-chain verification of transcripts, when instead of aggregating transcripts,
    /// valid payloads are aggregated.
    // pub fn verify_with_fs<D: EvaluationDomain<C::ScalarField>>(&self, params: &Params<C, D>, ss: SecretSharingWithWitness<C>) -> Result<DkgResult<C>, ()> {
    //     let mut fs = ark_transcript::Transcript::new_labeled(b"adkg_vrf::dkg::verifier");
    //     // TODO: hash the pks
    //     // fs.append(params);
    //     fs.append(&transcript);
    //     let rng = &mut fs.challenge(b"rng");
    //     if self.verify(params, &transcript, rng) {
    //         Ok(transcript.payload)
    //     } else {
    //         Err(())
    //     }
    // }

    // TODO: check params
    #[must_use]
    pub fn verify<R: Rng>(&self, ss: &SecretSharingWithWitness<C>, signer_pks: &[C::G2Affine], rng: &mut R) -> Result<(), ()> {
        let payload = &ss.payload;

        // 1, 2, 3, 4
        // Merges the equations from `Ceremony::verify_transcript_unoptimized` with random coefficients `r1, r2, r3`.
        let (r1, z) = (C::ScalarField::rand(rng), C::ScalarField::rand(rng));
        let r2 = r1.square();
        let r3 = r2 * r1;

        let _t = start_timer!(|| "Interpolation");
        let lis_size_n_at_z = self.domain_size_n.lagrange_basis_at(z);
        let (lis_size_t_at_z, lis_size_t_at_0) = {
            let mut lis_size_t_at_z = self.domain_size_t.lagrange_basis_at(z);
            let mut lis_size_t_at_0 = self.domain_size_t.lagrange_basis_at(C::ScalarField::zero());
            lis_size_t_at_z.resize(self.config.n, C::ScalarField::zero());
            lis_size_t_at_0.resize(self.config.n, C::ScalarField::zero());
            (lis_size_t_at_z, lis_size_t_at_0)
        };
        end_timer!(_t);

        let a_coeffs: Vec<_> = lis_size_n_at_z.iter()
            .zip(lis_size_t_at_z)
            .zip(lis_size_t_at_0)
            .map(|((li_n_z, li_t_z), li_t_0)| {
                (C::ScalarField::one() - r1) * li_n_z + r1 * li_t_z - r2 * li_t_0
            })
            .collect();

        let _t = start_timer!(|| "1xG1 + 2xG2 MSMs");
        let a_term = C::G1::msm(&ss.a, &a_coeffs)
            .map_or_else(|_err| Err(()), |x| Ok(x))?;
        let bgpk_at_z = C::G2::msm(&payload.bgpk, &lis_size_n_at_z)
            .map_or_else(|_err| Err(()), |x| Ok(x))?;
        let pk_at_z = C::G2::msm(&signer_pks, &lis_size_n_at_z)
            .map_or_else(|_err| Err(()), |x| Ok(x))?;
        end_timer!(_t);

        if C::multi_pairing(
            &[a_term + payload.c * r2 + payload.h1 * r3, -self.config.g1, payload.h1.into()],
            &[self.config.g2, bgpk_at_z + payload.h2 * r3, pk_at_z],
        ).is_zero() {
            Ok(())
        } else {
            Err(())
        }
    }
}

#[cfg(test)]
impl<C: Pairing> crate::pvss::Params<C> {
    pub fn verify_transcript_unoptimized<R: Rng>(&self, ss: &SecretSharingWithWitness<C>, rng: &mut R) {
        // 2. h2 has the same dlog as h1
        assert_eq!(C::pairing(ss.payload.h1, self.config.g2), C::pairing(self.config.g1, ss.payload.h2));
        // 3. `A`s are the evaluations of a degree `t` polynomial in the exponent
        self.verify_as(&ss, rng);
        // 4. `C = f(0).g1`
        self.verify_c(&ss);
        // 5. `bgpk`s are well-formed
        self.verify_bgpks(&ss, rng);
    }

    // Checks that `bgpk_j = f_i(w^j).g2 + sh_i.pk_j, j = 0,...,n-1`.
    // For that we interpolate 3 degree `< n` polynomials in the exponent:
    // 1. `bgpk(w^j).g2 = bgpk_j`,
    // 2. `f(w^j).g1 = A_j`, and
    // 3. `pk(w^j).g2 = pk_j`.
    // Then `bgpk(z) = f(z) + sh.pk(z)`, and, as `h1 = sh_i.g1`,
    // we can check that `e(g1, bgpk(z)) = e(f(z), g2) + e(h1, pk(z))`.
    fn verify_bgpks<R: Rng>(&self, ss: &SecretSharingWithWitness<C>, rng: &mut R) {
        let z = C::ScalarField::rand(rng);
        let lis_at_z = BarycentricDomain::of_size(self.config.domain, self.config.n).lagrange_basis_at(z);
        let f_at_z_g1 = C::G1::msm(&ss.a, &lis_at_z).unwrap();
        let bgpk_at_z_g2 = C::G2::msm(&ss.payload.bgpk, &lis_at_z).unwrap();
        let pk_at_z_g2 = C::G2::msm(&self.signer_pks, &lis_at_z).unwrap();
        let lhs = C::pairing(self.config.g1, bgpk_at_z_g2);
        let rhs = C::pairing(f_at_z_g1, self.config.g2) + C::pairing(ss.payload.h1, pk_at_z_g2);
        assert_eq!(lhs, rhs);
    }

    fn verify_as<R: Rng>(&self, t: &SecretSharingWithWitness<C>, rng: &mut R) {
        let z = C::ScalarField::rand(rng);
        let ls_deg_n_at_z = BarycentricDomain::of_size(self.config.domain, self.config.n).lagrange_basis_at(z);
        let ls_deg_t_at_z = BarycentricDomain::of_size(self.config.domain, self.config.t).lagrange_basis_at(z);
        let f_deg_n_at_z = C::G1::msm(&t.a, &ls_deg_n_at_z);
        let f_deg_t_at_z = C::G1::msm(&t.a[..self.config.t], &ls_deg_t_at_z);
        assert_eq!(f_deg_n_at_z, f_deg_t_at_z);
    }

    fn verify_c(&self, ss: &SecretSharingWithWitness<C>) {
        use ark_ec::CurveGroup;
        let ls_at_0 = BarycentricDomain::of_size(self.config.domain, self.config.t).lagrange_basis_at(C::ScalarField::zero());
        let f_at_0 = C::G1::msm(&ss.a[..self.config.t], &ls_at_0).unwrap();
        assert_eq!(ss.payload.c, f_at_0.into_affine());
    }
}
