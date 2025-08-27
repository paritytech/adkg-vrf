use ark_ff::{FftField, Field};
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, EvaluationDomain};
use ark_std::{end_timer, iter, start_timer};
use ark_std::{vec, vec::Vec};

/// Utilities for Lagrange interpolation over a subset of an FFT domain.

/// The vanishing polynomial of `x`.
/// `v(X) = X - x`
fn v<F: Field>(x: F) -> DensePolynomial<F> {
    DensePolynomial::from_coefficients_vec(vec![-x, F::one()]) // -x*0 + 1*X
}

/// The vanishing polynomial of the set `x_1,...,x_n`.
/// `v(X) = (X - x_1) * ... * (X - x_n)`
/// Runs in `O(nlog^2(n))` when the field `F` contains a radix-2 domain.
/// `n` doesn't have to be a power of `2`.
fn v_of<F: FftField>(xs: &[F]) -> DensePolynomial<F> {
    let n = xs.len();
    let mut product_tree = Vec::with_capacity(2 * n - 1);
    product_tree.extend(xs.iter().map(|&xi| v(xi)));
    for i in 0..n - 1 {
        product_tree.push(&product_tree[2 * i] * &product_tree[2 * i + 1])
    }
    product_tree[2 * n - 2].clone()
}

/// `f'`, the formal derivative of `f`.
fn diff<F: Field>(f: &DensePolynomial<F>) -> DensePolynomial<F> {
    let df_coeffs = f.iter()
        .enumerate()
        .skip(1)
        .map(|(i, ci)| F::from(i as u32) * ci)
        .collect();
    DensePolynomial::from_coefficients_vec(df_coeffs)
}

/// A set of interpolation points together with the precomputed weights.
/// The points don't have to form a multiplicative group.
/// After the weights are computed, evaluating an interpolant at a point takes `O(n)`.
#[derive(Clone)]
pub struct BarycentricDomain<F: Field> {
    /// Interpolation points `x_1,...,x_n`.
    xs: Vec<F>,
    /// Inverted barycentric weights of the points,
    /// `1 / w_j = prod(x_j - x_k, k != j), j = 1,...,n`.
    ws_inv: Vec<F>,
}

impl<F: FftField> BarycentricDomain<F> {

    /// The set of interpolation points is a subset of the `fft_domain` identified by the `bitmask`.
    /// The weights are computed as the evaluations of the derivative of the vanishing polynomial
    /// of the interpolation points over the FFT domain in `O(nlog^2(n))`.
    pub fn from_subset<D: EvaluationDomain<F>>(fft_domain: D, bitmask: &[bool]) -> Self {
        assert_eq!(bitmask.len(), fft_domain.size());
        let xs = fft_domain.elements()
            .zip(bitmask)
            .filter_map(|(xi, bi)| bi.then_some(xi))
            .collect::<Vec<_>>();
        let _t = start_timer!(|| "Product tree");
        let v = v_of(&xs);
        end_timer!(_t);
        let dv = diff(&v);
        let dv_over_domain = dv.evaluate_over_domain(fft_domain);
        let ws_inv = dv_over_domain.evals.into_iter()
            .zip(bitmask)
            .filter_map(|(xi, bi)| bi.then_some(xi))
            .collect::<Vec<_>>();
        Self {
            xs,
            ws_inv,
        }
    }

    /// The interpolation points are `x_k = w^k, k = 0,...,n-1`,
    /// where `w` is the generator of the `fft_domain`.
    pub fn of_size<D: EvaluationDomain<F>>(fft_domain: D, n: usize) -> Self {
        assert!(n <= fft_domain.size());
        let bitmask = {
            let mut bitmask = vec![true; n];
            bitmask.resize(fft_domain.size(), false);
            bitmask
        };
        Self::from_subset(fft_domain, &bitmask)
    }
}

impl<F: Field> BarycentricDomain<F> {

    /// Computes the weights of the set `x_1, ..., x_n` in `O(n^2)` using the formula
    /// `1 / w_j = prod(x_j - x_k, k != j), j = 1,...,n`.
    pub fn from_set(xs: Vec<F>) -> Self {
        let ws_inv = xs.iter().map(|xj| xs.iter()
            .filter_map(|xk| (xk != xj).then_some(*xj - xk))
            .product(),
        ).collect();
        Self {
            xs,
            ws_inv,
        }
    }

    /// Evaluates the Lagrange basis polynomials over the set `x_1, ..., x_n` at `z`.
    /// Returns `L_j(z), j = 1,...,n`, where `deg(L_j) = n - 1`, `L_j(x_j) = 1` and `L_j(x_k) = 0` for `k != j`.
    // Follows the formulae from https://people.maths.ox.ac.uk/trefethen/barycentric.pdf, p.3.
    // Fails if `x_j = x_k` for `j != k`.
    pub fn lagrange_basis_at(&self, z: F) -> Vec<F> {
        let (cs, l_at_z) = self._lagrange_basis_at(z);
        // L_j(z) = c_j.l(z)
        cs.iter()
            .map(|cj| l_at_z * cj)
            .collect()
    }

    pub fn evaluate(&self, ys: &[F], z: F) -> F {
        let (cs, l_at_z) = self._lagrange_basis_at(z);
        let p_at_z = l_at_z * cs.into_iter()
            .zip(ys)
            .map(|(cj, yj)| cj * yj)
            .sum::<F>();
        p_at_z
    }

    pub fn _lagrange_basis_at(&self, z: F) -> (Vec<F>, F) {
        // z - x_1, ..., z - x_n
        let z_minus_xs: Vec<F> = iter::repeat(z).zip(&self.xs)
            .map(|(z, xj)| z - xj)
            .collect();

        // l(z) = (z - x_1) ... (z - x_n)
        let l_at_z: F = z_minus_xs.iter()
            .product();

        // c_j = w_j / (z - x_j) = 1 / [(1 / w_j) * (z - x_j)]
        let cs = {
            let mut cs_inv: Vec<F> = self.ws_inv.iter().zip(z_minus_xs)
                .map(|(wj_inv, z_minus_xj)| z_minus_xj * wj_inv)
                .collect();
            ark_ff::batch_inversion(&mut cs_inv);
            cs_inv
        };

        (cs, l_at_z)
    }
}

pub fn powers<F: Field>(base: F) -> impl Iterator<Item=F> {
    iter::successors(Some(F::one()), move |prev| Some(base * prev))
}

#[cfg(test)]
mod tests {
    use ark_ec::CurveGroup;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::{DenseUVPolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial};
    use ark_std::rand::Rng;
    use ark_std::UniformRand;
    use ark_std::{end_timer, format, start_timer, test_rng};

    use super::*;

    #[test]
    // x_j = u^j, where u is a primary root of unity of order n
    fn over_subgroup() {
        let rng = &mut test_rng();

        let n = 16;
        let fft_domain = GeneralEvaluationDomain::new(n).unwrap();
        let d = BarycentricDomain::from_subset(fft_domain, &vec![true; n]);

        let z = ark_bls12_381::Fr::rand(rng);

        let ls1 = d.lagrange_basis_at(z);
        let ls2 = fft_domain.evaluate_all_lagrange_coefficients(z);
        assert_eq!(ls1, ls2);
    }

    #[test]
    // p(z) = p(x_1).L_1(z) + ... + p(x_n).L_n(z), if deg(p) < n
    fn over_random_set() {
        let rng = &mut test_rng();

        let n = 7;
        let xs = (0..n)
            .map(|_| ark_bls12_381::Fr::rand(rng))
            .collect::<Vec<_>>();
        let d = BarycentricDomain::from_set(xs.clone());

        let p = DensePolynomial::rand(n - 1, rng);
        let ys = xs.iter()
            .map(|xj|p.evaluate(xj))
            .collect::<Vec<_>>();
        let z = ark_bls12_381::Fr::rand(rng);

        let p_at_z = d.evaluate(&ys, z);

        assert_eq!(p_at_z, p.evaluate(&z));
    }

    #[test]
    fn test_vanishing_poly() {
        let n = 2usize.pow(10);
        let domain = GeneralEvaluationDomain::<ark_bls12_381::Fr>::new(n).unwrap();
        let xs = domain.elements().collect::<Vec<_>>();

        let _t = start_timer!(|| "Product tree");
        let v = v_of(&xs);
        end_timer!(_t);

        assert_eq!(v, domain.vanishing_polynomial().into());
    }

    #[test]
    // the weights computed with different formulas are the same
    fn test_barycentric_weights() {
        let rng = &mut test_rng();

        let (log_n, t) = (9, 2.0 / 3.0);
        let n = 2usize.pow(log_n);
        let domain = GeneralEvaluationDomain::<ark_bls12_381::Fr>::new(n).unwrap();
        let bitmask = _random_bits(n, t, rng);

        let _t = start_timer!(|| format!("Inverted barycentric weights, log(n)={}, t~{}", log_n, t));
        let d1 = BarycentricDomain::from_subset(domain, &bitmask);
        end_timer!(_t);

        let xs = d1.xs;
        let _t = start_timer!(|| format!("Naive inverted barycentric weights, log(n)={}, t~{}", log_n, t));
        let d2 = BarycentricDomain::from_set(xs);
        end_timer!(_t);

        assert_eq!(d1.ws_inv, d2.ws_inv);
    }

    fn _bench_barycentric_weights<F: FftField>(log_n: u32) {
        let n = 2usize.pow(log_n);
        let domain = GeneralEvaluationDomain::<F>::new(n).unwrap();
        let bitmask = vec![true; n];
        let _t = start_timer!(|| format!("Inverted barycentric weights, log(n)={}", log_n));
        let _d = BarycentricDomain::from_subset(domain, &bitmask);
        end_timer!(_t);
    }

    #[test]
    #[ignore]
    fn bench_barycentric_weights() {
        _bench_barycentric_weights::<ark_bls12_381::Fr>(10);
        _bench_barycentric_weights::<ark_bls12_381::Fr>(16);
        _bench_barycentric_weights::<ark_bls12_381::Fr>(20);
    }

    fn _bench_msm<G: CurveGroup>(log_n: u32) {
        let rng = &mut test_rng();

        let n = 2usize.pow(log_n);
        let bases = (0..n).map(|_| G::rand(rng)).collect::<Vec<_>>();
        let bases_affine = G::normalize_batch(&bases);
        let exps = (0..n).map(|_| G::ScalarField::rand(rng)).collect::<Vec<_>>();

        let _t = start_timer!(|| format!("MSM, log(n)={}", log_n));
        let _msm = G::msm(&bases_affine, &exps);
        end_timer!(_t);
    }

    #[test]
    #[ignore]
    fn bench_against_msm() {
        _bench_barycentric_weights::<ark_bls12_381::Fr>(10);
        _bench_msm::<ark_bls12_381::G1Projective>(10);
        _bench_msm::<ark_bls12_381::G2Projective>(10);
    }

    fn _random_bits<R: Rng>(n: usize, density: f64, rng: &mut R) -> Vec<bool> {
        (0..n).map(|_| rng.gen_bool(density)).collect()
    }
}
