use ark_ec::pairing::{Pairing, PairingOutput};
use ark_ec::PrimeGroup;

use crate::bls::vanilla::StandaloneSig;
use crate::pvss::SecretSharing;

// Used to verify aggregated threshold signatures.
// `c = f(0).g1` is the public key associated with the dealing.
// `(h1, h2)` bind signers to the dealing via `bgpk_j = gsk_j + sk_j.h2`.
pub struct ThresholdVk<C: Pairing> {
    c: C::G1Affine,
    h1: C::G1Affine,
    h2: C::G2Affine,
    // todo: skip serialization
    g1: C::G1,
    g2: C::G2,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AggThresholdSig<C: Pairing> {
    pub(crate) bls_sig_with_pk: StandaloneSig<C>,
    pub(crate) bgpk: C::G2Affine,
}

impl<C: Pairing> ThresholdVk<C> {
    pub fn from_share(share: &SecretSharing<C>) -> Self {
        //TODO: consume?
        Self {
            c: share.c,
            h1: share.h1,
            h2: share.h2,
            // TODO:
            g1: C::G1::generator(),
            g2: C::G2::generator(),
        }
    }

    pub fn verify_unoptimized(&self, sig: &AggThresholdSig<C>, message: C::G1) {
        sig.bls_sig_with_pk
            .verify_unoptimized(message, self.g2.into());
        assert_eq!(
            C::pairing(self.g1.into(), sig.bgpk),
            C::multi_pairing(
                &[self.c, self.h1],
                &[self.g2.into(), sig.bls_sig_with_pk.pk]
            )
        );
    }

    pub fn vuf_unoptimized(&self, sig: &AggThresholdSig<C>, message: C::G1) -> PairingOutput<C> {
        self.verify_unoptimized(sig, message);
        C::multi_pairing(
            &[(-message).into(), sig.bls_sig_with_pk.sig],
            &[sig.bgpk, self.h2],
        )
    }
}
