use ark_ec::pairing::Pairing;
use ark_ec::CurveGroup;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::vec::Vec;

use crate::koe;
use crate::old_dkg::DkgResult;

/// Standalone or aggregated transcript with the witness.
// TODO: add weights?
#[derive(Clone, Debug, PartialEq, Eq, CanonicalDeserialize, CanonicalSerialize)]
pub struct DkgTranscript<C: Pairing> {
    pub payload: DkgResult<C>,

    // witness data
    /// Commitment to the secret polynomial `A_j = f(w^j).g1, j = 0,...,n-1`
    pub(crate) a: Vec<C::G1Affine>,
    /// Proofs of knowledge of the exponents `(f_i(0), sh_i)`
    /// such that `C_i=f_i(0).g1` and `h1_i=sh_i.g1` for every dealer `i = 1,...,k`.
    pub(crate) koe_proofs: Vec<KoeProof<C>>,
}

/// Proof that the dealer `i` knows her secrets.
#[derive(Clone, Debug, PartialEq, Eq, CanonicalDeserialize, CanonicalSerialize)]
pub(crate) struct KoeProof<C: Pairing> {
    /// `C_i = f_i(0).g1`
    pub(crate) c_i: C::G1Affine,
    /// `h1_i = sh_i.g1`
    pub(crate) h1_i: C::G1Affine,
    /// `s_i` is a proof of knowledge of the discrete logs of `(C_i, h1_i)` with respect to `g1`.
    pub(crate) koe_proof: koe::Proof<C::G1>,
}

impl<C: Pairing> DkgTranscript<C> {
    pub fn merge_with(self, others: &[Self]) -> Self {
        let mut others = others.to_vec();
        others.push(self);
        Self::merge(&others)
    }

    pub fn merge(transcripts: &[Self]) -> Self {
        let n = transcripts[0].a.len();
        let a = (0..n).map(|j| {
            transcripts.iter()
                .map(|t| t.a[j])
                .sum::<C::G1>().into_affine()
        }).collect();

        let payload = transcripts.iter()
            .map(|t| t.payload.clone())
            .collect::<Vec<_>>();
        let payload = DkgResult::merge(&payload);

        let koe_proofs = transcripts.iter()
            .flat_map(|t| t.koe_proofs.clone())
            .collect::<Vec<_>>();

        Self {
            payload,
            a,
            koe_proofs,
        }
    }
}
