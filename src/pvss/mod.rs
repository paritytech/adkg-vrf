mod dealer;
pub mod verifier;

pub use verifier::Verifier;

use ark_ec::pairing::Pairing;
use ark_ec::{CurveGroup, PrimeGroup};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};

/// An Aggregatable **Publicly Verifiable Secret Sharing** (PVSS) protocol.

/// The purpose of a *threshold secret sharing* protocol is to decompose a secret into `n` chunks
/// in such a way that any `t` (but not less than `t`) of them is enough to reconstruct the secret.
/// A PVSS additionally requires the consistency of all the chunks to be publicly verifiable.
/// That dictates to output all the chunks at once, each encrypted to its recipient.
/// Finally, in an *aggregatable PVSS*, `2` valid outputs `pvss(a)` and `pvss(b)`
/// can be combined into a valid output for the sum of inputs `pvss(a + b)`.

/// In our case the secret is an elliptic curve point. The recipients are BLS signers,
/// whose public keys `(pk_1,...,pk_n)` in `G2` are known in order. `1 <= t <= n`.
/// TODO:

#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct Config<
    C: Pairing,
    D: EvaluationDomain<C::ScalarField> = GeneralEvaluationDomain<<C as Pairing>::ScalarField>,
> {
    /// The number of signers.
    pub n: usize,
    /// The threshold, i.e. the minimal number of signers required to reconstruct the shared secret.
    pub t: usize,
    /// An FFT-friendly multiplicative subgroup of the field, of size not less than `n`.
    /// The evaluation points are the first `n` elements of the subgroup: `x_j = w^j, j = 0,...,n-1`,
    /// where `w` is the generator of the subgroup.
    pub domain: D,
    /// Generator of G1.
    pub g1: C::G1,
    /// Generator of G2.
    pub g2: C::G2,
}

impl<C: Pairing> Config<C> {
    pub fn new(n: usize, t: usize) -> Result<Self, ()> {
        if !(n > 0 && t > 0 && t <= n) {
            // todo: test t = 1, t = n
            return Err(());
        }
        let domain = GeneralEvaluationDomain::new(n).ok_or(())?;
        Ok(Self {
            n,
            t,
            domain,
            g1: C::G1::generator(),
            g2: C::G2::generator(),
        })
    }
}

/// Parameters of a PVSS.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Params<
    C: Pairing,
    D: EvaluationDomain<C::ScalarField> = GeneralEvaluationDomain<<C as Pairing>::ScalarField>,
> {
    pub config: Config<C, D>,
    /// The signers' bls public keys in G2.
    /// **Proofs of possession must be checked for these keys.**
    pub signer_pks: Vec<C::G2Affine>,
}

impl<C: Pairing> Params<C> {
    pub fn new(signer_pks: Vec<C::G2Affine>, t: usize) -> Result<Self, ()> {
        let n = signer_pks.len();
        let config = Config::new(n, t)?;
        Ok(Self { config, signer_pks })
    }
}

/// A dealer samples a scalar `sh` and a degree `t-1` polynomial `f` with scalar coefficients.
/// 1. `sh` contributes a pair of points `(h1 = sh.g1, h2 = sh.g2)` in `G1xG2` with equal DLOGs.
/// 2. `f(0).g2` is the secret to share among the signers. The corresponding public key `c = f(0).g1`.
///    `(f(w^j).g2, j = 0,...,n-1)` are the shares of the secret. With any`t` of them one can interpolate `f(0).g2`.
///     Each share `f(w^j).g2` is encrypted to the corresponding signer `pk_j` with ElGamal as
///    `(bgpk_j = f(w^j).g2 + sh.pk_j, j = 0,...,n-1)`. A signer `j` can (though doesn't need to)
///     use the BLS secret key `sk_j` to decrypt `(bgpk_j, h2)` as `f(w^j).g2 = bgpk_j - sk_j.h2`.

/// Shares of a secret, and the public keys.
#[derive(Clone, Debug, PartialEq, Eq, CanonicalDeserialize, CanonicalSerialize)]
pub struct SecretSharing<C: Pairing> {
    /// The public key corresponding to the shared secret key.
    /// `c = f(0).g1`
    pub c: C::G1Affine,
    /// Shares of the secret, encrypted to the signers.
    /// `(bgpk_j = f(w^j).g2 + sh.pk_j, j = 0,...,n-1)`
    pub bgpk: Vec<C::G2Affine>,
    /// `h1 = sh.g1`
    pub h1: C::G1Affine,
    /// `h2 = sh.g2`
    pub h2: C::G2Affine,
}

/// `SecretSharing` with the corresponding witness attached.
#[derive(Clone, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct SecretSharingWithWitness<C: Pairing> {
    /// Witness: commitment to the secret polynomial `f`.
    /// `(a_j = f(w^j).g1, j = 0,...,n-1)`
    pub a: Vec<C::G1Affine>,
    pub payload: SecretSharing<C>,
}

impl<C: Pairing> SecretSharingWithWitness<C> {
    pub fn aggregate_with(self, mut others: Vec<Self>) -> Self {
        others.push(self);
        Self::aggregate(&others)
    }

    pub fn aggregate(sharings: &[Self]) -> Self {
        let n = sharings[0].a.len();
        let agg_ss = SecretSharing {
            c: sharings
                .iter()
                .map(|s| s.payload.c)
                .sum::<C::G1>()
                .into_affine(),
            h1: sharings
                .iter()
                .map(|s| s.payload.h1)
                .sum::<C::G1>()
                .into_affine(),
            h2: sharings
                .iter()
                .map(|s| s.payload.h2)
                .sum::<C::G2>()
                .into_affine(),
            bgpk: (0..n)
                .map(|j| {
                    sharings
                        .iter()
                        .map(|s| s.payload.bgpk[j])
                        .sum::<C::G2>()
                        .into_affine()
                })
                .collect(),
        };
        SecretSharingWithWitness {
            a: (0..n)
                .map(|j| sharings.iter().map(|s| s.a[j]).sum::<C::G1>().into_affine())
                .collect(),
            payload: agg_ss,
        }
    }
}
