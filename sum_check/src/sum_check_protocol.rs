use ark_bn254::Fq;
use ark_ff::{BigInteger, PrimeField};
use fiat_shamir::fiat_shamir_transcript::Transcript;
use multilinear_polynomial::multilinear_polynomial_evaluation::MultilinearPoly;

#[derive(Debug)]
struct Proof {
    proof_polynomials: Vec<Vec<Fq>>,
    claimed_sum: Fq,
}

fn get_round_partial_polynomial_proof(polynomial: &[Fq]) -> Vec<Fq> {
    let mid_point = polynomial.len() / 2;
    let (zeros, ones) = polynomial.split_at(mid_point);

    let poly_proof = vec![zeros.iter().sum(), ones.iter().sum()];

    poly_proof
}

fn prove(polynomial: &MultilinearPoly<Fq>) -> Proof {
    let mut transcript = Transcript::<Fq>::new();
    transcript.append(&fq_vec_to_bytes(&polynomial.evaluation));

    let claimed_sum: Fq = polynomial.evaluation.iter().sum();
    transcript.append(&fq_vec_to_bytes(&[claimed_sum]));

    let num_rounds = polynomial.evaluation.len().ilog2();
    let mut proof_polynomials = Vec::with_capacity(num_rounds as usize);
    let mut current_poly = polynomial.clone();

    for _ in 0..num_rounds {
        let proof_poly = get_round_partial_polynomial_proof(&current_poly.evaluation);

        transcript.append(&fq_vec_to_bytes(&proof_poly));

        proof_polynomials.push(proof_poly);

        let random_challenge = transcript.get_random_challenge();

        current_poly = current_poly.partial_evaluate(0, &random_challenge);
    }

    Proof {
        proof_polynomials,
        claimed_sum,
    }
}

fn verify(polynomial: &MultilinearPoly<Fq>, proof: Proof) -> bool {
    let mut transcript = Transcript::<Fq>::new();
    transcript.append(&fq_vec_to_bytes(&polynomial.evaluation));
    transcript.append(&fq_vec_to_bytes(&[proof.claimed_sum]));

    let mut current_poly = polynomial.clone();
    let mut random_challenges = Vec::with_capacity(proof.proof_polynomials.len());
    let mut expected_sum = proof.claimed_sum;

    for poly in proof.proof_polynomials {
        let poly = MultilinearPoly::new(poly.to_vec());

        if poly.evaluation.iter().sum::<Fq>() != expected_sum {
            return false;
        }

        transcript.append(&fq_vec_to_bytes(&poly.evaluation));
        let random_challenge = transcript.get_random_challenge();

        expected_sum =
            poly.evaluation[0] + random_challenge * (poly.evaluation[1] - poly.evaluation[0]);

        current_poly = current_poly.partial_evaluate(0, &random_challenge);

        random_challenges.push(random_challenge);
    }

    let poly_eval_sum = polynomial.evaluate(random_challenges);
    expected_sum == poly_eval_sum
}

fn fq_vec_to_bytes(values: &[Fq]) -> Vec<u8> {
    values
        .iter()
        .flat_map(|x| x.into_bigint().to_bytes_le())
        .collect()
}

#[cfg(test)]
mod test {
    use ark_bn254::Fq;
    use multilinear_polynomial::multilinear_polynomial_evaluation::MultilinearPoly;

    use crate::sum_check_protocol::{prove, verify, Proof};

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

        let proof = prove(&initial_polynomial);

        let is_verified = verify(&initial_polynomial, proof);

        assert_eq!(is_verified, true);
    }

    #[test]
    fn test_invalid_proof_doesnt_verify() {
        let initial_polynomial =
            MultilinearPoly::new(vec![Fq::from(0), Fq::from(3), Fq::from(2), Fq::from(5)]);

        let false_proof = Proof {
            claimed_sum: Fq::from(20),
            proof_polynomials: vec![
                vec![Fq::from(3), Fq::from(9)],
                vec![Fq::from(1), Fq::from(2)],
            ],
        };

        let is_verified = verify(&initial_polynomial, false_proof);

        assert_eq!(is_verified, false);
    }
}
