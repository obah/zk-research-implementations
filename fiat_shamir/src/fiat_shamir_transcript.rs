use ark_ff::PrimeField;
use sha3::{Digest, Keccak256};
use std::marker::PhantomData;

pub struct Transcript<F: PrimeField> {
    _field: PhantomData<F>,
    hasher: Keccak256,
}

impl<F: PrimeField> Transcript<F> {
    pub fn new() -> Self {
        Self {
            _field: PhantomData,
            hasher: Keccak256::new(),
        }
    }

    pub fn append(&mut self, preimage: &[u8]) {
        self.hasher.update(preimage)
    }

    pub fn get_random_challenge(&self) -> F {
        let final_hash = self.hasher.clone().finalize(); //? this might be problematic

        F::from_le_bytes_mod_order(&final_hash)
    }
}

#[cfg(test)]
mod test {
    use super::Transcript;
    use ark_bn254::Fq;

    #[test]
    fn it_hashes() {
        let mut transcript: Transcript<Fq> = Transcript::new();
        transcript.append("zero knowledge".as_bytes());

        let random_challenge = transcript.get_random_challenge();

        dbg!(random_challenge);
    }
}
