use ark_ff::PrimeField;

use crate::multilinear_polynomial_evaluation::MultilinearPoly;

pub struct ProductPoly<F: PrimeField> {
    evaluation: Vec<MultilinearPoly<F>>,
}

impl<F: PrimeField> ProductPoly<F> {
    pub fn new(evaluations: Vec<Vec<F>>) -> Self {
        let polys = evaluations
            .iter()
            .map(|evaluation| MultilinearPoly::new(evaluation.to_vec()))
            .collect();

        Self { evaluation: polys }
    }

    fn evaluate(&self, values: Vec<F>) -> F {
        self.evaluation
            .iter()
            .map(|poly| poly.evaluate(values.clone()))
            .product()
    }

    pub fn partial_evaluate(&self, value: &F) -> Vec<Vec<F>> {
        self.evaluation
            .iter()
            .map(|poly| {
                let partial_res = poly.partial_evaluate(0, value);

                partial_res.evaluation
            })
            .collect()
    }
}

#[cfg(test)]
mod test {
    use ark_bn254::Fq;

    use super::ProductPoly;

    #[test]
    fn it_evaluates_multiple_polys() {
        let evaluations = vec![
            vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(3)],
            vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)],
        ];

        let product_polys = ProductPoly::new(evaluations);

        let values = vec![Fq::from(2), Fq::from(3)];

        let expected_evaluation = Fq::from(216);

        let result = product_polys.evaluate(values);

        assert_eq!(expected_evaluation, result);
    }

    #[test]
    fn it_partially_evaluates_multiple_polys() {
        let evaluations = vec![
            vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(3)],
            vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)],
        ];

        let product_polys = ProductPoly::new(evaluations);

        let value = Fq::from(2);

        let expected_evaluation = vec![
            vec![Fq::from(0), Fq::from(6)],
            vec![Fq::from(0), Fq::from(4)],
        ];

        let result = product_polys.partial_evaluate(&value);

        assert_eq!(result, expected_evaluation);
    }
}
