use ark_ec::pairing::Pairing;
use ark_ec::{CurveGroup, ScalarMul};
use ark_ff::Zero;
use ark_poly::{DenseUVPolynomial, Polynomial};

use crate::pvss::{Params, SecretSharing, SecretSharingWithWitness};
use ark_poly::univariate::DensePolynomial;
use ark_std::rand::Rng;
use ark_std::{end_timer, start_timer, UniformRand};

impl<C: Pairing> Params<C> {
    pub fn deal<R: Rng>(&self, rng: &mut R) -> SecretSharingWithWitness<C> {
        let f = DensePolynomial::rand(self.config.t - 1, rng);
        let sh = C::ScalarField::rand(rng);
        self._deal(f, sh)
    }

    /// Samples a random degree `t-1` polynomial `f` with constant term `ssk` (in other words, `f(0) = ssk`)
    /// and shares the secret `ssk.g2` using `f`. `(h1 = sh.g1, h2 = sh.g2)`.
    pub fn deal_secrets<R: Rng>(&self, ssk: C::ScalarField, sh: C::ScalarField, rng: &mut R) -> SecretSharingWithWitness<C> {
        let t = self.config.t;
        let mut coeffs = Vec::with_capacity(t);
        coeffs.push(ssk); // constant term
        coeffs.extend(&DensePolynomial::rand(t - 2, rng).coeffs); // ensures the leading coeff is not `0`
        assert_eq!(coeffs[0], ssk);
        assert!(!coeffs[t - 1].is_zero());
        let f = DensePolynomial::from_coefficients_vec(coeffs);
        self._deal(f, sh)
    }

    fn _deal(&self, f_mon: DensePolynomial<C::ScalarField>, sh: C::ScalarField) -> SecretSharingWithWitness<C> {
        let ssk = f_mon[0];
        assert!(!ssk.is_zero());
        assert!(!sh.is_zero());
        assert_eq!(f_mon.degree(), self.config.t - 1);
        let f_lag: Vec<C::ScalarField> = f_mon.evaluate_over_domain(self.config.domain)
            .evals.into_iter()
            .take(self.config.n)
            .collect();

        // For log_n = 10, storing the precomputed tables would save only 5% of the total dealing time.
        let _t = start_timer!(|| "Commitment to the secret polynomial in G1 and G2");
        let f_lag_g1 = self.config.g1.batch_mul(&f_lag); // `f(w^j).g1, j = 0,...,n-1`, the pvss witness `a`
        let f_lag_g2 = self.config.g2.batch_mul(&f_lag); // `f(w^j).g2, j = 0,...,n-1`, shares of the secret
        end_timer!(_t);

        let _t = start_timer!(|| "Key exchange");
        let shared_secrets: Vec<_> = self.signer_pks.iter()
            .map(|&pk_j| pk_j * sh)
            .collect();
        end_timer!(_t);

        // And now we just add the shared secrets to make the secret shares secret.
        let secret_shares = f_lag_g2;
        let bgpk: Vec<_> = secret_shares.into_iter()
            .zip(shared_secrets)
            .map(|(share, delta)| share + delta)
            .collect();
        let bgpk = C::G2::normalize_batch(&bgpk);

        // Can be batched, but who cares.
        let c = (self.config.g1 * ssk).into_affine();
        let h1 = (self.config.g1 * sh).into_affine();
        let h2 = (self.config.g2 * sh).into_affine();

        SecretSharingWithWitness {
            a: f_lag_g1,
            payload: SecretSharing { c, bgpk, h1, h2 },
        }
    }
}