use ark_ec::pairing::Pairing;
use ark_ec::{CurveGroup, ScalarMul};
use ark_poly::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_poly::EvaluationDomain;
use ark_std::rand::Rng;
use ark_std::{end_timer, start_timer, UniformRand};
use ark_std::{vec, vec::Vec};

use crate::koe;
use crate::old_dkg::transcript::{DkgTranscript, KoeProof};
use crate::old_dkg::{Ceremony, DkgResult};

impl<'a, C: Pairing, D: EvaluationDomain<C::ScalarField>> Ceremony<'a, C, D> {
    //TODO: cryptorng?
    pub fn deal<R: Rng>(&self, rng: &mut R) -> DkgTranscript<C> {
        // dealer's secrets
        let (f_mon, sh) = (DensePolynomial::rand(self.t - 1, rng), C::ScalarField::rand(rng));
        let ssk = f_mon[0];
        let f_lag: Vec<C::ScalarField> = f_mon.evaluate_over_domain(self.domain).evals.into_iter()
            .take(self.n)
            .collect();

        // For log_n = 10, storing the precomputed tables would save only 5% of the total dealing time.
        let _t = start_timer!(|| "Commitment to the secret polynomial in G1 and G2");
        let f_lag_g1 = self.g1.batch_mul(&f_lag);
        let f_lag_g2 = self.g2.batch_mul(&f_lag);
        end_timer!(_t);

        let _t = start_timer!(|| "Key exchange");
        let shared_keys: Vec<_> = self.bls_pks.iter()
            .map(|&pk_j| pk_j * sh)
            .collect();
        end_timer!(_t);

        let bgpk: Vec<_> = f_lag_g2.into_iter()
            .zip(shared_keys)
            .map(|(f_j, sk_j)| f_j + sk_j)
            .collect();
        let bgpk = C::G2::normalize_batch(&bgpk);

        let c = self.g1 * ssk;
        let h1 = self.g1 * sh;
        let h2 = self.g2 * sh;

        let instance = koe::Instance { base: self.g1, points: vec![c, h1] };
        let statement = koe::Statement { instance, dlogs: vec![ssk, sh] };
        let koe_proof = statement.prove(rng);

        // Can be batched, but who cares.
        let c = c.into_affine();
        let h1 = h1.into_affine();
        let h2 = h2.into_affine();
        let payload = DkgResult { bgpk, h1, h2, c };
        let koe_proof = KoeProof { c_i: c, h1_i: h1, koe_proof };
        DkgTranscript {
            payload,
            a: f_lag_g1,
            koe_proofs: vec![koe_proof],
        }
    }
}

#[cfg(test)]
mod tests {
    use ark_poly::GeneralEvaluationDomain;
    use ark_std::{end_timer, format, start_timer, test_rng};

    use super::*;

    fn _bench_dkg<C: Pairing>(f: usize) {
        let rng = &mut test_rng();

        // n = 3f+1 -- number of signers,
        // t = 2f+1 -- threshold number of signers
        // k = f+1 -- number of dealers
        let (n, t, k) = (3 * f + 1, 2 * f + 1, f + 1);
        let signers = (0..n)
            .map(|_| C::G2Affine::rand(rng))
            .collect::<Vec<_>>();
        let params = Ceremony::<C, GeneralEvaluationDomain<C::ScalarField>>::setup(t, &signers);
        let _t = start_timer!(|| format!("Transcript generation, n = {}, t = {}", n, t));
        let transcript = params.deal(rng);
        end_timer!(_t);

        let _t = start_timer!(|| format!("Precomputation for transcript validation"));
        let transcript_verifier = params.verifier();
        end_timer!(_t);

        let _t = start_timer!(|| format!("Standalone transcript validation"));
        assert!(transcript_verifier.verify(&params, &transcript, rng));
        end_timer!(_t);

        let transcripts = vec![transcript; k];
        let _t = start_timer!(|| format!("Merging {} transcripts", k));
        let agg_transcript = DkgTranscript::merge(&transcripts);
        end_timer!(_t);

        let _t = start_timer!(|| format!("Aggregate transcript validation, k = {}", k));
        assert!(transcript_verifier.verify(&params, &agg_transcript, rng));
        end_timer!(_t);
    }

    #[cfg(all(not(debug_assertions), feature = "print-trace"))]
    fn _bench_batch_mul<C: CurveGroup>(log_n: u32) {
        let rng = &mut test_rng();

        let g = C::generator();
        let scalars: Vec<_> = (0..2usize.pow(log_n))
            .map(|_| C::ScalarField::rand(rng))
            .collect();

        let _t = start_timer!(|| format!("Batch multiplication, log(n) = {}", log_n));

        let _t_table = start_timer!(|| "Precomputation");
        let table = ark_ec::scalar_mul::BatchMulPreprocessing::new(g, scalars.len());
        end_timer!(_t_table);

        ark_std::println!("Table size = {}x{} affine points", table.table.len(), table.table[0].len());

        let _t_mul = start_timer!(|| "Multiplication");
        C::batch_mul_with_preprocessing(&table, &scalars);
        end_timer!(_t_mul);

        end_timer!(_t);
    }

    #[test]
    #[ignore]
    #[cfg(all(not(debug_assertions), feature = "print-trace"))]
    fn bench_batch_mul() {
        ark_std::println!("G1:");
        _bench_batch_mul::<ark_bls12_381::G1Projective>(10);
        ark_std::println!("G2:");
        _bench_batch_mul::<ark_bls12_381::G2Projective>(10);
    }

    #[test]
    #[ignore]
    #[cfg(all(not(debug_assertions), feature = "print-trace"))]
    // cargo test bench_dkg_10 --release --features="print-trace" -- --nocapture --ignored
    fn bench_dkg_10() {
        assert_eq!((2usize.pow(10) - 1) / 3, 341);
        _bench_dkg::<ark_bls12_381::Bls12_381>(341);
    }

    #[test]
    #[ignore]
    #[cfg(all(not(debug_assertions), feature = "print-trace"))]
    fn bench_dkg_16() {
        assert_eq!((2usize.pow(16) - 1) / 3, 21845);
        _bench_dkg::<ark_bls12_381::Bls12_381>(21845);
    }
}