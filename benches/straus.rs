use ark_ec::scalar_mul::glv::GLVConfig;
use ark_ec::short_weierstrass::Affine;
use ark_ec::CurveGroup;
use ark_std::{test_rng, UniformRand};
use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};

use adkg_vrf::straus;

fn short_msm<C: CurveGroup>(c: &mut Criterion) {
    let rng = &mut test_rng();
    let n = 4; // msm size
    let m = 100; // pre-generated input size

    let full_scalars = (0..m)
        .map(|_| C::ScalarField::rand(rng))
        .collect::<Vec<_>>();
    let small_scalars = (0..m)
        .map(|_| C::ScalarField::from(u128::rand(rng)))
        .collect::<Vec<_>>();
    let points = (0..m).map(|_| C::rand(rng)).collect::<Vec<_>>();
    let points = C::normalize_batch(&points);

    let mut i = 0;
    c.bench_with_input(BenchmarkId::new("Straus full scalar MSM", n), &n, |b, n| {
        b.iter_batched(
            || {
                i = (i + 1) % (m - n);
                (&points[i..i + n], &full_scalars[i..i + n])
            },
            |(points, scalars)| straus::short_msm(points, scalars),
            BatchSize::SmallInput,
        )
    });

    let mut i = 0;
    c.bench_with_input(
        BenchmarkId::new("Straus 128-bit scalar MSM", n),
        &n,
        |b, n| {
            b.iter_batched(
                || {
                    i = (i + 1) % (m - n);
                    (&points[i..i + n], &small_scalars[i..i + n])
                },
                |(points, scalars)| straus::short_msm(points, scalars),
                BatchSize::SmallInput,
            )
        },
    );
}

fn glv_vs_straus<C: GLVConfig>(c: &mut Criterion) {
    let rng = &mut test_rng();
    let mut bg = c.benchmark_group("Small MSMs");

    let m = 100; // pre-generated input size
    let full_scalars = (0..m)
        .map(|_| C::ScalarField::rand(rng))
        .collect::<Vec<_>>();
    let points = (0..m).map(|_| Affine::<C>::rand(rng)).collect::<Vec<_>>();

    for n in [2, 4] {
        let mut i = 0;
        bg.bench_with_input(BenchmarkId::new("Straus full scalar MSM", n), &n, |b, n| {
            b.iter_batched(
                || {
                    i = (i + 1) % (m - n);
                    (&points[i..i + n], &full_scalars[i..i + n])
                },
                |(points, scalars)| straus::short_msm(points, scalars),
                BatchSize::SmallInput,
            )
        });

        let mut i = 0;
        bg.bench_with_input(BenchmarkId::new("GLV-Straus MSM", n), &n, |b, n| {
            b.iter_batched(
                || {
                    i = (i + 1) % (m - n);
                    (&points[i..i + n], &full_scalars[i..i + n])
                },
                |(points, scalars)| straus::short_msm_glv(points, scalars),
                BatchSize::SmallInput,
            )
        });
    }
}

criterion_group!(
    benches,
    short_msm::<ark_bls12_381::G1Projective>,
    glv_vs_straus::<ark_bls12_381::g1::Config>,
);
criterion_main!(benches);
