use std::marker::PhantomData;

use ark_ff::PrimeField;
use ark_std::rand::rngs::StdRng;
use dusk_bls12_381::BlsScalar;
use dusk_poseidon::{Domain, Hash};
use rand::{Rng, SeedableRng};

pub struct Transcript<'a, F> {
    _field: PhantomData<F>,
    hasher: Hash<'a>,
}

impl<'a, F: PrimeField> Transcript<'a, F> {
    pub fn new() -> Self {
        Self {
            _field: PhantomData,
            hasher: Hash::new(Domain::Other),
        }
    }

    fn append(&mut self, preimage: &'a [BlsScalar]) {
        self.hasher.update(&preimage)
    }

    pub fn get_random_challenge(&mut self) -> F {
        let final_hash = self.hasher.finalize();

        let hash_bytes = final_hash
            .iter()
            .flat_map(|scalar| scalar.to_bytes())
            .collect::<Vec<u8>>();

        F::from_le_bytes_mod_order(&hash_bytes)
    }
}

fn bytes_to_scalars(bytes: &[u8]) -> Result<Vec<BlsScalar>, &'static str> {
    const SCALAR_SIZE: usize = 32;

    if bytes.len() % SCALAR_SIZE != 0 {
        return Err("Input byte slice length must be a multiple of Scalar size (32 bytes)");
    }

    let scalars = bytes
        .chunks(SCALAR_SIZE)
        .map(|chunk| {
            let array: [u8; SCALAR_SIZE] = chunk
                .try_into()
                .expect("Chunk size is guaranteed to be SCALAR_SIZE");
            BlsScalar::from_bytes(&array).expect("Invalid Scalar bytes")
        })
        .collect();

    Ok(scalars)
}
