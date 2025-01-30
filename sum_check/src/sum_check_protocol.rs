use ark_bn254::Fq;
use ark_ff::{BigInteger, PrimeField};
use fiat_shamir::fiat_shamir_transcript::Transcript;
use multilinear_polynomial::multilinear_polynomial_evaluation::MultilinearPoly;

struct Proof {
    initial_poly: MultilinearPoly<Fq>,
    claimed_sum: Fq,
    proof_polys: Vec<Vec<Fq>>,
    sum_proofs: Vec<Fq>,
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
            sum_proofs: Vec::new(),
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

    fn get_proof(&mut self) -> Proof {
        let sum_proof = self.get_sum_proof();
        let proof_poly = self.get_partial_polynomial_proof();

        let mut new_sum_proofs = self.proof.sum_proofs.clone();
        new_sum_proofs.push(sum_proof);

        let mut new_proof_polys = self.proof.proof_polys.clone();
        new_proof_polys.push(proof_poly.evaluation);

        Proof {
            initial_poly: self.polynomial.clone(),
            claimed_sum: sum_proof,
            proof_polys: new_proof_polys,
            sum_proofs: new_sum_proofs,
        }
    }

    fn verify_proof(&mut self, proof: Proof) -> bool {
        if proof.proof_polys.last().iter().len() == 1 {
            let mut random_challenges = Vec::new();

            for (i, poly) in proof.proof_polys.iter().enumerate() {
                let poly_bytes = fq_vec_to_bytes(poly);
                let sum_bytes = fq_vec_to_bytes(&vec![proof.sum_proofs[i]]);

                self.transcript.append(&poly_bytes);
                self.transcript.append(&sum_bytes);

                let random_challenge = self.transcript.get_random_challenge();

                random_challenges.push(random_challenge);

                let multi_poly = MultilinearPoly::new(poly.to_vec());

                multi_poly.partial_evaluate(random_challenge, 0);
            }

            let poly_eval_sum = self.polynomial.evaluate(random_challenges);

            proof.proof_polys.last().unwrap()[0] == poly_eval_sum
        } else {
            let last_poly = MultilinearPoly::new(proof.proof_polys.last().unwrap().to_vec());
            let eval_0 = last_poly.evaluate(vec![Fq::from(0)]);
            let eval_1 = last_poly.evaluate(vec![Fq::from(1)]);

            let poly_sum = eval_0 + eval_1;

            *proof.sum_proofs.last().unwrap() == poly_sum
        }
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
