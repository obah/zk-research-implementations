use crate::gkr_circuit::Operation;
use ark_ff::PrimeField;

fn add_mul_polynomials<F: PrimeField>(poly_a: Vec<F>, poly_b: Vec<F>, op: Operation) -> Vec<F> {
    let new_eval_len = poly_a.len() * poly_b.len();
    let mut new_eval = Vec::with_capacity(new_eval_len);

    for a in poly_a {
        for b in &poly_b {
            new_eval.push(op.apply(a, *b))
        }
    }

    new_eval
}

#[cfg(test)]
mod test {
    use super::add_mul_polynomials;
    use crate::gkr_circuit::Operation;
    use ark_bn254::Fq;

    #[test]
    fn it_add_polys_correctly() {
        let poly_a = vec![Fq::from(0), Fq::from(2)];
        let poly_b = vec![Fq::from(0), Fq::from(3)];

        let expected_poly = vec![Fq::from(0), Fq::from(3), Fq::from(2), Fq::from(5)];

        let result = add_mul_polynomials(poly_a, poly_b, Operation::Add);

        assert_eq!(result, expected_poly);

        let poly_a = vec![Fq::from(0), Fq::from(3)];
        let poly_b = vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)];

        let expected_poly = vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(2),
            Fq::from(3),
            Fq::from(3),
            Fq::from(3),
            Fq::from(5),
        ];

        let result = add_mul_polynomials(poly_a, poly_b, Operation::Add);

        assert_eq!(result, expected_poly);

        // let poly_a = vec![Fq::from(0), Fq::from(3), Fq::from(2), Fq::from(5)];
        // let poly_b = vec![Fq::from(0), Fq::from(4), Fq::from(1), Fq::from(5)];

        // let expected_poly = vec![
        //     Fq::from(0),
        //     Fq::from(3),
        //     Fq::from(4),
        //     Fq::from(7),
        //     Fq::from(1),
        //     Fq::from(4),
        //     Fq::from(5),
        //     Fq::from(8),
        //     Fq::from(2),
        //     Fq::from(5),
        //     Fq::from(6),
        //     Fq::from(9),
        //     Fq::from(3),
        //     Fq::from(6),
        //     Fq::from(7),
        //     Fq::from(10),
        // ];

        // let result = add_polys(poly_a, poly_b);

        // assert_eq!(result, expected_poly);
    }

    #[test]
    fn it_multiplies_polys_correctly() {
        let poly_a = vec![Fq::from(0), Fq::from(2)];
        let poly_b = vec![Fq::from(0), Fq::from(3)];

        let expected_poly = vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(6)];

        let result = add_mul_polynomials(poly_a, poly_b, Operation::Mul);

        assert_eq!(result, expected_poly);

        let poly_a = vec![Fq::from(0), Fq::from(3)];
        let poly_b = vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)];

        let expected_poly = vec![
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(6),
        ];

        let result = add_mul_polynomials(poly_a, poly_b, Operation::Mul);

        assert_eq!(result, expected_poly);
    }
}
