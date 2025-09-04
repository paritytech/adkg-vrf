use ark_ec::scalar_mul::glv::GLVConfig;
use ark_ec::short_weierstrass::{Affine, Projective};
use ark_ec::{AffineRepr, CurveGroup};
use ark_ff::{AdditiveGroup, BigInteger, BitIteratorBE, PrimeField, Zero};
use ark_std::{end_timer, format, start_timer};
use ark_std::{iter, vec, vec::Vec};

/// Simplest form of Straus' multi-scalar multiplication.
/// See the Handbook of Elliptic and Hyperelliptic Curve Cryptography,
/// [Algorithm 9.23](https://hyperelliptic.org/HEHCC/chapters/chap09.pdf).

/// Precomputes sums of all subsets of the list of `points`.
// The ordering is as in `Self::bits_to_index`.
fn table<C: AffineRepr>(points: &[C], w: u32) -> Vec<C> {
    let c = 2usize.pow(w);
    let mut table = vec![C::Group::zero()];
    for p in points {
        // P, 2P, ..., (c-1)P, where c = 2^w
        let multiples_of_p: Vec<C::Group> =
            iter::successors(Some(p.into_group()), move |prev| Some(*p + *prev))
                .take(c - 1)
                .collect();
        // TODO: batchconvert to affine?
        let new_rows: Vec<C::Group> = multiples_of_p
            .iter()
            .flat_map(|&kp| table.iter().map(move |&prev_row| prev_row + kp))
            .collect();
        table.extend(new_rows)
    }
    C::Group::normalize_batch(&table)
}

fn bits_to_digit<I: Iterator<Item = bool>>(bits: I, powers_of_2: &[u32]) -> u32 {
    bits.zip(powers_of_2.iter())
        .filter_map(|(bit, power)| bit.then_some(power))
        .sum::<u32>()
}

fn digits_to_index<I: Iterator<Item = u32>>(digits: I, powers_of_c: &[u32]) -> usize {
    digits
        .zip(powers_of_c.iter())
        .map(|(digit, power)| digit * power)
        .sum::<u32>() as usize
}

/// Converts `bits` highlighting a subset of the points to the index at which the sum of the subset is located in the table.
// The powers of `2` should start from `1`, so the least significant bit goes first.
// fn bits_to_index<I: Iterator<Item=bool>>(bits: I, powers_of_2: &[u32]) -> usize {
//     bits.zip(powers_of_2.iter())
//         .filter_map(|(bit, power)| bit.then_some(power))
//         .sum::<u32>() as usize
// }

/// Pads the binary decomposition of the `scalar` by `0`s from the left (msbf),
/// so that the bit length is divided evenly by the window size `w`.
fn to_msbf_bits_padded<F: PrimeField>(scalar: F, w: usize) -> Vec<bool> {
    let repr_bit_len = F::BigInt::NUM_LIMBS * 64;
    let extra_bits = repr_bit_len % w;
    let padding_len = if extra_bits == 0 { 0 } else { w - extra_bits };
    iter::repeat(false)
        .take(padding_len)
        .chain(BitIteratorBE::new(scalar.into_bigint()))
        .collect()
}

fn to_base_c_digits<F: PrimeField>(scalars: &[F], w: usize) -> Vec<Vec<u32>> {
    let powers_of_2 = iter::successors(Some(1u32), move |prev| Some(prev << 1))
        .take(w)
        .collect::<Vec<_>>();

    scalars
        .iter()
        .map(|&s| {
            to_msbf_bits_padded(s, w)
                .chunks(w)
                .map(|w_bit_chunk| bits_to_digit(w_bit_chunk.iter().rev().cloned(), &powers_of_2))
                .collect()
        })
        .collect()
}

/// Converts a list of `scalars` into a list of lookup indices per window of size `w`.
fn indices<F: PrimeField>(scalars: &[F], w: usize) -> Vec<usize> {
    // let c = 2usize.pow(w);
    // represent each scalar as a list of base `c` digits
    let scalars_base_c = to_base_c_digits(scalars, w);
    let powers_of_c = iter::successors(Some(1u32), move |prev| Some(prev << w))
        .take(scalars.len())
        .collect::<Vec<_>>();
    // let repr_bit_len = F::BigInt::NUM_LIMBS * 64;
    let skip = 0;
    // let scalar_bit_len = F::MODULUS_BIT_SIZE as usize;
    // let skip = repr_bit_len - scalar_bit_len;
    (skip..scalars_base_c[0].len())
        .map(|i| {
            let slice = scalars_base_c.iter().map(|s| s[i]);
            digits_to_index(slice, &powers_of_c)
        })
        .collect()
}

pub fn short_msm<C: AffineRepr>(points: &[C], scalars: &[C::ScalarField]) -> C::Group {
    short_msm_windowed(points, scalars, 1)
}

pub fn short_msm_windowed<C: AffineRepr>(
    points: &[C],
    scalars: &[C::ScalarField],
    w: usize,
) -> C::Group {
    let _t_table = start_timer!(|| "Points");
    let table = table(points, w as u32);
    end_timer!(_t_table);
    let _t_indices = start_timer!(|| "Indices");
    let indices = indices(scalars, w);
    end_timer!(_t_indices);
    let mut acc = C::Group::zero();
    for idx in indices.into_iter().skip_while(|&idx| idx == 0) {
        for _ in 0..w {
            acc.double_in_place();
        }
        acc += table[idx]
    }
    acc
}

fn glv_decomposition<C: GLVConfig>(
    p: Affine<C>,
    s: C::ScalarField,
) -> Vec<(Affine<C>, C::ScalarField)> {
    let ((sgn_s1, s1), (sgn_s2, s2)) = C::scalar_decomposition(s);
    let mut p2 = C::endomorphism_affine(&p);
    let mut p1 = p;
    if !sgn_s1 {
        p1 = -p1;
    }
    if !sgn_s2 {
        p2 = -p2;
    }
    vec![(p1, s1), (p2, s2)]
}

pub fn short_msm_glv<C: GLVConfig>(
    points: &[Affine<C>],
    scalars: &[C::ScalarField],
) -> Projective<C> {
    assert_eq!(points.len(), scalars.len());
    let _t_glv = start_timer!(|| format!("{} GLV decompositions", points.len()));
    let (glv_points, glv_scalars): (Vec<Affine<C>>, Vec<C::ScalarField>) = points
        .iter()
        .zip(scalars)
        .flat_map(|(p, s)| glv_decomposition(*p, *s))
        .unzip();
    end_timer!(_t_glv);
    let _t_msm = start_timer!(|| format!("straus {}-msm", glv_points.len()));
    let res = short_msm(&glv_points, &glv_scalars);
    end_timer!(_t_msm);
    res
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_std::println;
    use ark_std::{test_rng, UniformRand};

    fn _test_short_msm<C: GLVConfig>(n: usize) {
        let rng = &mut test_rng();

        let scalars = (0..n)
            .map(|_| C::ScalarField::rand(rng))
            .collect::<Vec<_>>();
        let points = (0..n).map(|_| Affine::<C>::rand(rng)).collect::<Vec<_>>();

        let _t_naive = start_timer!(|| format!("Naive {}-msm", n));
        let res: Projective<C> = points.iter().zip(scalars.iter()).map(|(&p, s)| p * s).sum();
        end_timer!(_t_naive);
        println!();

        let _t_straus_1 = start_timer!(|| format!("Straus msm, n = {}, w = 1", n));
        let res1 = short_msm(&points, &scalars);
        end_timer!(_t_straus_1);
        assert_eq!(res1, res);
        println!();

        let _t_straus_2 = start_timer!(|| format!("Straus msm, n = {}, w = 2", n));
        let res2 = short_msm_windowed(&points, &scalars, 2);
        end_timer!(_t_straus_2);
        assert_eq!(res2, res);
        println!();

        let _t_straus_3 = start_timer!(|| format!("Straus msm, n = {}, w = 3", n));
        let res3 = short_msm_windowed(&points, &scalars, 3);
        end_timer!(_t_straus_3);
        assert_eq!(res3, res);
        println!();

        let _t_glv = start_timer!(|| format!("GLV {}-msm", n));
        let res4 = short_msm_glv(&points, &scalars);
        end_timer!(_t_glv);
        assert_eq!(res4, res);
    }

    #[test]
    fn test_short_msm() {
        _test_short_msm::<ark_bls12_381::g1::Config>(2);
        _test_short_msm::<ark_bls12_381::g1::Config>(3);
        _test_short_msm::<ark_bls12_381::g1::Config>(4);
    }
}
