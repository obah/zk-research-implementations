use ark_bn254::Fq;
use ark_ff::{BigInteger, PrimeField};
use fiat_shamir::fiat_shamir_transcript::Transcript;
use multilinear_polynomial::multilinear_polynomial_evaluation::MultilinearPoly;

struct Proof {
    initial_poly: MultilinearPoly<Fq>,
    claimed_sum: Fq,
    proof_polys: Vec<Vec<Fq>>,
}

struct SumCheck {
    polynomial: MultilinearPoly<Fq>,
    transcript: Transcript<Fq>,
    proof: Proof,
}

//// ! check and confirm the polynomial in sumcheck isnt being mutated anywhere

impl SumCheck {
    fn new(poly: MultilinearPoly<Fq>) -> Self {
        let mut new_transcript = Transcript::new();

        let poly_bytes = fq_vec_to_bytes(&poly.evaluation);

        new_transcript.append(&poly_bytes);

        let new_proof = Proof {
            initial_poly: poly.clone(),
            claimed_sum: Fq::from(0),
            proof_polys: Vec::new(),
        };

        Self {
            polynomial: poly,
            transcript: new_transcript,
            proof: new_proof,
        }
    }

    fn get_sum_proof(&mut self) -> Fq {
        let sum_proof: Fq = self.polynomial.evaluation.iter().sum();

        let sum_bytes = fq_vec_to_bytes(&vec![sum_proof]);

        self.transcript.append(&sum_bytes);

        sum_proof
    }

    fn get_partial_polynomial_proof(&self) -> MultilinearPoly<Fq> {
        let mid_point = self.polynomial.evaluation.len() / 2;

        let (zeros, ones) = self.polynomial.evaluation.split_at(mid_point);

        MultilinearPoly::new(vec![zeros.iter().sum(), ones.iter().sum()])
    }

    fn get_proof(&self) -> Proof {
        let sum_proof = self.get_sum_proof();
        let proof_poly = self.get_partial_polynomial_proof();

        let new_proof_polys = self.proof.proof_polys;
        new_proof_polys.push(proof_poly);

        // Proof {
        //     initial_poly: self.polynomial,
        //     claimed_sum: ,
        //     proof_polys: new_proof_polys,
        // }

        todo!()
    }

    fn verify_sum_proof(poly_proof: Self, sum_proof: Fq) -> bool {
        //rework this this to include oracle check

        let eval_0 = poly_proof.polynomial.evaluate(vec![Fq::from(0)]);
        let eval_1 = poly_proof.polynomial.evaluate(vec![Fq::from(1)]);

        let poly_sum = eval_0 + eval_1;

        sum_proof == poly_sum
    }

    fn initiate_next_round(&self) -> MultilinearPoly<Fq> {
        let random_challenge = self.transcript.get_random_challenge();

        let new_poly = self.polynomial.partial_evaluate(random_challenge, 0);

        new_poly
    }
}

fn fq_vec_to_bytes(value: &Vec<Fq>) -> Vec<u8> {
    let value_bytes: Vec<u8> = value
        .iter()
        .flat_map(|x| x.into_bigint().to_bytes_le())
        .collect();

    value_bytes
}

// Proof {initial_poly, claimed_sum, uni_poly}
// Prover::new(poly)
// Prover::prove(initial_poly) -> Proof
// Verfier::new()
// Verifier::verify(Proof, initial_poly)
// Transcript::new()
