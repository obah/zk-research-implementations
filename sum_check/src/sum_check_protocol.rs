use ark_bn254::Fq;
use ark_ff::{BigInteger, PrimeField};
use fiat_shamir::fiat_shamir_transcript::Transcript;
use multilinear_polynomial::multilinear_polynomial_evaluation::MultilinearPoly;

struct Proof {
    proof_polynomials: Vec<Vec<Fq>>,
    claimed_sums: Vec<Fq>,
}

struct Prover {
    transcript: Transcript<Fq>,
}

struct Verifier {
    transcript: Transcript<Fq>,
}

impl Prover {
    fn new(polynomial: &MultilinearPoly<Fq>) -> Self {
        let mut new_transcript = Transcript::new();

        let polynomial_bytes = fq_vec_to_bytes(&polynomial.evaluation);

        new_transcript.append(&polynomial_bytes);

        Self {
            transcript: new_transcript,
        }
    }

    fn get_sum_proof(&mut self, polynomial: &Vec<Fq>) -> Fq {
        let sum_proof: Fq = polynomial.iter().sum();

        let sum_bytes = fq_vec_to_bytes(&vec![sum_proof]);

        self.transcript.append(&sum_bytes);

        sum_proof
    }

    fn get_partial_polynomial_proof(&mut self, polynomial: &Vec<Fq>) -> Vec<Fq> {
        let mid_point = polynomial.len() / 2;

        let (zeros, ones) = polynomial.split_at(mid_point);

        let poly_proof = vec![zeros.iter().sum(), ones.iter().sum()];

        let poly_proof_bytes = fq_vec_to_bytes(&poly_proof);

        self.transcript.append(&poly_proof_bytes);

        poly_proof
    }

    fn proof(&mut self, polynomial: &MultilinearPoly<Fq>) -> Proof {
        let num_of_rounds = polynomial.evaluation.len().ilog2();

        let mut claimed_sums = Vec::<Fq>::new();
        let mut proof_polynomials = Vec::<Vec<Fq>>::new();
        let mut current_poly = polynomial.clone();

        for _ in 0..num_of_rounds {
            //evaluate first one and
            let claimed_sum = self.get_sum_proof(&current_poly.evaluation);
            claimed_sums.push(claimed_sum);

            let proof_poly = self.get_partial_polynomial_proof(&current_poly.evaluation);
            proof_polynomials.push(proof_poly);

            let random_challenge = self.transcript.get_random_challenge();

            current_poly = current_poly.partial_evaluate(random_challenge, 0);
        }

        Proof {
            proof_polynomials,
            claimed_sums,
        }
    }
}

impl Verifier {
    fn new(polynomial: &MultilinearPoly<Fq>) -> Self {
        let mut new_transcript = Transcript::new();

        let polynomial_bytes = fq_vec_to_bytes(&polynomial.evaluation);

        new_transcript.append(&polynomial_bytes);

        Self {
            transcript: new_transcript,
        }
    }

    fn verify(&mut self, polynomial: &MultilinearPoly<Fq>, proof: Proof) -> bool {
        let mut random_challenges = Vec::new();

        for (i, poly) in proof.proof_polynomials.iter().enumerate() {
            let poly_bytes = fq_vec_to_bytes(poly);
            let sum_bytes = fq_vec_to_bytes(&vec![proof.claimed_sums[i]]);

            self.transcript.append(&poly_bytes);
            self.transcript.append(&sum_bytes);

            let random_challenge = self.transcript.get_random_challenge();

            random_challenges.push(random_challenge);

            let multi_poly = MultilinearPoly::new(poly.to_vec());

            multi_poly.partial_evaluate(random_challenge, 0);
        }

        let poly_eval_sum = polynomial.evaluate(random_challenges);

        proof.proof_polynomials.last().unwrap()[0] == poly_eval_sum
    }
}

fn fq_vec_to_bytes(value: &Vec<Fq>) -> Vec<u8> {
    let value_bytes: Vec<u8> = value
        .iter()
        .flat_map(|x| x.into_bigint().to_bytes_le())
        .collect();

    value_bytes
}

#[cfg(test)]
mod test {
    use super::{Proof, Prover, Verifier};
    use ark_bn254::Fq;
    use multilinear_polynomial::multilinear_polynomial_evaluation::MultilinearPoly;

    #[test]
    fn test_valid_proving_and_verification() {
        let initial_polynomial = MultilinearPoly::new(vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(0),
            Fq::from(10),
            Fq::from(0),
            Fq::from(17),
        ]);

        let mut prover = Prover::new(&initial_polynomial);
        let proof = prover.proof(&initial_polynomial);

        let mut verifier = Verifier::new(&initial_polynomial);
        let verified = verifier.verify(&initial_polynomial, proof);

        assert_eq!(verified, true);
    }

    #[test]
    fn test_invalid_proof_doesnt_verify() {
        let initial_polynomial =
            MultilinearPoly::new(vec![Fq::from(0), Fq::from(3), Fq::from(2), Fq::from(5)]);

        let false_proof = Proof {
            claimed_sums: vec![Fq::from(20), Fq::from(9)],
            proof_polynomials: vec![
                vec![Fq::from(3), Fq::from(9)],
                vec![Fq::from(1), Fq::from(2)],
            ],
        };

        let mut verifier = Verifier::new(&initial_polynomial);
        let verified = verifier.verify(&initial_polynomial, false_proof);

        assert_eq!(verified, false);
    }
}
