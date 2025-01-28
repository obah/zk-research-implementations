use crate::implementations::multilinear_polynomial::MultilinearPoly;
use ark_bn254::Fq;

use super::fiat_shamir::Transcipt;

struct SumCheck {
    polynomial: MultilinearPoly<Fq>,
    transcript: Transcipt<Fq>,
}

impl SumCheck {
    fn new(poly: MultilinearPoly<Fq>) -> Self {
        let new_transcript = Transcipt::new(&poly.evaluation);

        SumCheck {
            polynomial: poly,
            transcript: new_transcript,
        }
    }

    fn get_sum_proof(self) -> Fq {
        self.polynomial.evaluation.iter().sum()
    }

    fn get_partial_polynomial_proof(self) -> Self {
        let mut result = self.polynomial;

        let bits = result.evaluation.len().ilog2() - 1;
        let values = vec![Fq::from(0), Fq::from(1)];

        for i in 1..bits {
            let mut temp_poly = MultilinearPoly::new(vec![Fq::from(0)]);

            for value in &values {
                let partial_result = result.partial_evaluate(*value, i.try_into().unwrap());

                temp_poly += partial_result;
            }

            result = temp_poly;
        }

        SumCheck::new(result)
    }

    fn verify_sum_proof(poly_proof: Self, sum_proof: Fq) -> bool {
        let eval_0 = poly_proof.polynomial.evaluate(vec![Fq::from(0)]);
        let eval_1 = poly_proof.polynomial.evaluate(vec![Fq::from(1)]);

        let poly_sum = eval_0 + eval_1;

        sum_proof == poly_sum
    }

    fn initiate_next_round(mut self) -> MultilinearPoly<Fq> {
        let random_challenge = self.transcript.get_random_challenge();

        let new_poly = self.polynomial.partial_evaluate(random_challenge, 0);

        self.transcript.add_polynomial(new_poly.evaluation.clone());

        new_poly
    }
}
