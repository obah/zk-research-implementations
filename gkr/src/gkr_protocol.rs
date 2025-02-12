use crate::gkr_circuit::Operation;

use ark_bn254::Fq;

use fiat_shamir::fiat_shamir_transcript::Transcript;
use multilinear_polynomial::{
    multilinear_poly_product::{ProductPoly, SumPoly},
    multilinear_polynomial_evaluation::MultilinearPoly,
};

fn add_mul_polynomials(poly_a: &[Fq], poly_b: &[Fq], op: Operation) -> MultilinearPoly<Fq> {
    let new_eval_len = poly_a.len() * poly_b.len();
    let mut new_eval = Vec::with_capacity(new_eval_len);

    for a in poly_a {
        for b in poly_b {
            new_eval.push(op.apply(*a, *b))
        }
    }

    MultilinearPoly::new(new_eval)
}

fn get_fbc_poly(
    random_challenge: Fq,
    add_i: MultilinearPoly<Fq>,
    mul_i: MultilinearPoly<Fq>,
    w_b: &[Fq],
    w_c: &[Fq],
) -> SumPoly<Fq> {
    let new_add_i = add_i.partial_evaluate(0, &random_challenge);
    let new_mul_i = mul_i.partial_evaluate(0, &random_challenge);

    let summed_w_poly = add_mul_polynomials(w_b, w_c, Operation::Add);
    let multiplied_w_poly = add_mul_polynomials(w_b, w_c, Operation::Mul);

    let add_w_eval = vec![new_add_i.evaluation, summed_w_poly.evaluation];
    let mul_w_eval = vec![new_mul_i.evaluation, multiplied_w_poly.evaluation];

    let add_eval_product = ProductPoly::new(add_w_eval);
    let mul_eval_product = ProductPoly::new(mul_w_eval);

    SumPoly::new(vec![add_eval_product, mul_eval_product])
}

fn prove() {
    let mut transcript = Transcript::<Fq>::new();
    //? if output is one, i.e array len = 0, push a 0 to it
    // prover sends output poly to verifier (transcript)
    // prover also sends claimed sum to verifier (transcript)
}

fn verify() {}

#[cfg(test)]
mod test {
    use super::{add_mul_polynomials, get_fbc_poly};
    use crate::gkr_circuit::Operation;
    use ark_bn254::Fq;
    use multilinear_polynomial::{
        multilinear_poly_product::{ProductPoly, SumPoly},
        multilinear_polynomial_evaluation::MultilinearPoly,
    };

    #[test]
    fn it_add_polys_correctly() {
        let poly_a = &[Fq::from(0), Fq::from(2)];
        let poly_b = &[Fq::from(0), Fq::from(3)];

        let expected_poly = vec![Fq::from(0), Fq::from(3), Fq::from(2), Fq::from(5)];

        let result = add_mul_polynomials(poly_a, poly_b, Operation::Add);

        assert_eq!(result.evaluation, expected_poly);

        let poly_a = &[Fq::from(0), Fq::from(3)];
        let poly_b = &[Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)];

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

        assert_eq!(result.evaluation, expected_poly);
    }

    #[test]
    fn it_multiplies_polys_correctly() {
        let poly_a = &[Fq::from(0), Fq::from(2)];
        let poly_b = &[Fq::from(0), Fq::from(3)];

        let expected_poly = vec![Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(6)];

        let result = add_mul_polynomials(poly_a, poly_b, Operation::Mul);

        assert_eq!(result.evaluation, expected_poly);

        let poly_a = &[Fq::from(0), Fq::from(3)];
        let poly_b = &[Fq::from(0), Fq::from(0), Fq::from(0), Fq::from(2)];

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

        assert_eq!(result.evaluation, expected_poly);
    }

    #[test]
    fn test_get_fbc_poly() {
        let r_c = Fq::from(5);
        let add_i = MultilinearPoly::new(vec![
            Fq::from(0),
            Fq::from(1),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
        ]);
        let mul_i = MultilinearPoly::new(vec![Fq::from(0); 8]);
        let w_1_poly = &[Fq::from(2), Fq::from(12)];

        let fbc_poly = get_fbc_poly(r_c, add_i.clone(), mul_i.clone(), w_1_poly, w_1_poly);

        let one = ProductPoly::new(vec![
            add_i.evaluation,
            vec![Fq::from(4), Fq::from(14), Fq::from(14), Fq::from(24)],
        ]);
        let two = ProductPoly::new(vec![
            mul_i.evaluation,
            vec![Fq::from(4), Fq::from(24), Fq::from(24), Fq::from(144)],
        ]);

        let expected_result = SumPoly::new(vec![one, two]);

        assert_eq!(fbc_poly.polys, expected_result.polys);
    }
}
