use ark_bls12_381::{Bls12_381, Fr, G1Projective as G1, G2Projective as G2};
use ark_ec::{pairing::Pairing, AffineRepr, PrimeGroup, ScalarMul};
use ark_ff::PrimeField;
use multilinear_polynomial::multilinear_polynomial_evaluation::{MultilinearPoly, Operation};

struct Proof {
    quotients: Vec<G1>,
}

struct KZG {
    polynomial: MultilinearPoly<Fr>,
    g_1: G1,
    g_2: G2,
    g2_taus: Vec<G2>,
    g1_lagrange_basis: Vec<G1>,
}

impl KZG {
    fn new(polynomial: &MultilinearPoly<Fr>, taus: Vec<Fr>) -> Self {
        if taus.len() != polynomial.num_of_vars {
            panic!("invalid taus or polynomials");
        }

        let g_1 = G1::generator();
        let g_2 = G2::generator();

        let (g2_taus, g1_lagrange_basis) = KZG::run_trusted_setup(polynomial, g_1, g_2, taus);

        Self {
            polynomial: polynomial.clone(),
            g_1,
            g_2,
            g2_taus,
            g1_lagrange_basis,
        }
    }

    fn run_trusted_setup(
        poly: &MultilinearPoly<Fr>,
        g_1: G1,
        g_2: G2,
        taus: Vec<Fr>,
    ) -> (Vec<G2>, Vec<G1>) {
        let g2_taus_affine = g_2.batch_mul(&taus);

        let g2_taus: Vec<G2> = g2_taus_affine
            .into_iter()
            .map(|point| point.into_group())
            .collect();

        let lagrange_basis = get_lagrange_basis(poly.num_of_vars, &taus, g_1);

        (g2_taus, lagrange_basis)
    }

    fn commit(&self) -> G1 {
        evaluate_poly_with_l_basis_in_g1(&self.polynomial.evaluation, &self.g1_lagrange_basis)
    }

    fn open(&self, opening_values: &[Fr]) -> Fr {
        self.polynomial.evaluate(opening_values.to_vec())
    }

    fn get_proof(&self, opened_value: Fr, opening_values: &[Fr]) -> Proof {
        let mut poly_minus_v = MultilinearPoly::new(
            self.polynomial
                .evaluation
                .iter()
                .map(|eval| *eval - opened_value)
                .collect::<Vec<_>>(),
        );

        let mut q_i: Vec<G1> = Vec::with_capacity(opening_values.len());
        let full_n_vars = self.polynomial.num_of_vars;
        for value in opening_values {
            let mut quotient = get_quotient(&poly_minus_v, 0);

            let mut quotient_vars = quotient.len().ilog2();

            while quotient_vars < (full_n_vars as u32) {
                quotient = blow_up_poly(&quotient, quotient.len() * 2);
                quotient_vars += 1;
            }

            let quotient_eval =
                evaluate_poly_with_l_basis_in_g1(&quotient, &self.g1_lagrange_basis);

            q_i.push(quotient_eval);

            let remainder_poly = get_remainder(&poly_minus_v, *value, 0);

            poly_minus_v = MultilinearPoly::new(remainder_poly);
        }

        Proof { quotients: q_i }
    }

    fn verify(
        &self,
        commitment: G1,
        opened_value: Fr,
        proof: Proof,
        opening_values: &[Fr],
    ) -> bool {
        if proof.quotients.len() != opening_values.len() {
            panic!("num of quotients in proof not equal to num of opening values");
        }

        let lhs = commitment - self.g_1.mul_bigint(opened_value.into_bigint());
        let lhs_gt = Bls12_381::pairing(lhs, self.g_2);

        let mut rhs = Vec::with_capacity(opening_values.len());

        for i in 0..opening_values.len() {
            let quotient = proof.quotients[i];
            let opening_value_g2 = self.g_2.mul_bigint(opening_values[i].into_bigint());
            let factor = self.g2_taus[i] - opening_value_g2;

            let rhs_i = Bls12_381::pairing(quotient, factor);
            rhs.push(rhs_i);
        }

        let rhs_gt = rhs.iter().sum();

        lhs_gt == rhs_gt
    }
}

fn evaluate_poly_with_l_basis_in_g1(poly_evaluations: &[Fr], lagrange_basis: &[G1]) -> G1 {
    if poly_evaluations.len() != lagrange_basis.len() {
        panic!("invalid polynomial or lagrange basis");
    }

    poly_evaluations
        .iter()
        .zip(lagrange_basis.iter())
        .map(|(a, b)| b.mul_bigint(a.into_bigint()))
        .sum()
}

fn get_remainder(poly: &MultilinearPoly<Fr>, value: Fr, bit: usize) -> Vec<Fr> {
    poly.partial_evaluate(bit, &value).evaluation
}

fn get_quotient(poly: &MultilinearPoly<Fr>, bit: usize) -> Vec<Fr> {
    let mut eval_0 = poly.partial_evaluate(bit, &Fr::from(0));
    let mut eval_1 = poly.partial_evaluate(bit, &Fr::from(1));

    if eval_0.evaluation.len() > eval_1.evaluation.len() {
        eval_1 = MultilinearPoly::new(blow_up_poly(&eval_1.evaluation, eval_0.evaluation.len()));
    } else if eval_0.evaluation.len() < eval_1.evaluation.len() {
        eval_0 = MultilinearPoly::new(blow_up_poly(&eval_0.evaluation, eval_1.evaluation.len()));
    }

    (eval_1 - eval_0).evaluation
}

fn blow_up_poly(poly: &[Fr], bigger_poly_len: usize) -> Vec<Fr> {
    let blow_up_factor = bigger_poly_len / poly.len();

    let blow_up_poly = vec![Fr::from(1); blow_up_factor];

    MultilinearPoly::tensor_add_mul_polynomials(&blow_up_poly, poly, Operation::Mul).evaluation
}

fn generate_bhc(bits: usize) -> Vec<Vec<u8>> {
    let total = 1usize << bits;

    (0..total)
        .map(|i| {
            (0..bits)
                .map(|j| ((i >> (bits - 1 - j)) & 1) as u8)
                .collect()
        })
        .collect()
}

fn get_lagrange_basis(num_of_vars: usize, unenc_taus: &[Fr], g_1: G1) -> Vec<G1> {
    if num_of_vars < 1 {
        panic!("Invalid num of vars for lagrange basis");
    }

    let mut result = Vec::with_capacity(1 << num_of_vars);
    let bhc = generate_bhc(num_of_vars);

    for layer in bhc {
        let mut layer_eval = Fr::from(1);

        for (i, bit) in layer.iter().enumerate() {
            let bit_scalar;
            if *bit == 0 {
                bit_scalar = Fr::from(1) - unenc_taus[i];
            } else {
                bit_scalar = unenc_taus[i];
            }

            layer_eval *= bit_scalar;
        }

        result.push(layer_eval);
    }

    g_1.batch_mul(&result)
        .into_iter()
        .map(|point| point.into_group())
        .collect()
}

#[cfg(test)]
mod test {
    use ark_bls12_381::{Fr, G1Projective};
    use ark_ec::{PrimeGroup, ScalarMul};
    use ark_ff::PrimeField;

    use super::*;

    #[test]
    fn test_blow_up_poly() {
        let poly = vec![Fr::from(0), Fr::from(4)];

        let expected_poly = vec![Fr::from(0), Fr::from(4), Fr::from(0), Fr::from(4)];

        let result = blow_up_poly(&poly, expected_poly.len());

        assert_eq!(expected_poly, result);
    }

    #[test]
    fn test_get_lagrange_basis() {
        let n_vars = 3;
        let unenc_taus = &[Fr::from(5), Fr::from(2), Fr::from(3)];
        let g_1 = G1Projective::generator();

        let lagrange_basis = get_lagrange_basis(n_vars, unenc_taus, g_1);

        let expected_l_basis = vec![
            Fr::from(-8),
            Fr::from(12),
            Fr::from(16),
            Fr::from(-24),
            Fr::from(10),
            Fr::from(-15),
            Fr::from(-20),
            Fr::from(30),
        ];

        let expected_l_basis_g1 = g_1.batch_mul(&expected_l_basis);

        assert_eq!(lagrange_basis, expected_l_basis_g1);
    }

    #[test]
    fn test_evaluate_poly_with_l_basis() {
        let poly_evals = &[
            Fr::from(0),
            Fr::from(4),
            Fr::from(0),
            Fr::from(4),
            Fr::from(0),
            Fr::from(4),
            Fr::from(3),
            Fr::from(7),
        ];

        let n_vars = 3;
        let unenc_taus = &[Fr::from(5), Fr::from(2), Fr::from(3)];
        let g_1 = G1Projective::generator();

        let lagrange_basis = get_lagrange_basis(n_vars, unenc_taus, g_1);

        let eval_result = evaluate_poly_with_l_basis_in_g1(poly_evals, &lagrange_basis);

        let expected_eval_result = g_1.mul_bigint(Fr::from(42).into_bigint());

        assert_eq!(eval_result, expected_eval_result);
    }

    #[test]
    fn test_get_remainder() {
        let main_poly = MultilinearPoly::new(vec![
            Fr::from(-72),
            Fr::from(-68),
            Fr::from(-54),
            Fr::from(-50),
        ]);
        let value = Fr::from(4);

        let remainder_result = get_remainder(&main_poly, value, 0);

        let expected_remainder = vec![Fr::from(0), Fr::from(4)];

        assert_eq!(remainder_result, expected_remainder);
    }

    #[test]
    fn test_get_quotient() {
        let main_poly = MultilinearPoly::new(vec![
            Fr::from(-72),
            Fr::from(-68),
            Fr::from(-54),
            Fr::from(-50),
        ]);

        let quotient_result = get_quotient(&main_poly, 0);

        let expected_quotient = Fr::from(18);

        assert_eq!(quotient_result[0], expected_quotient);
    }

    #[test]
    fn test_commit() {
        let poly_evals = &[
            Fr::from(0),
            Fr::from(4),
            Fr::from(0),
            Fr::from(4),
            Fr::from(0),
            Fr::from(4),
            Fr::from(3),
            Fr::from(7),
        ];
        let unenc_taus = vec![Fr::from(5), Fr::from(2), Fr::from(3)];
        let kzg_instance = KZG::new(&MultilinearPoly::new(poly_evals.to_vec()), unenc_taus);

        let commit_result = kzg_instance.commit();

        let expected_commit_result = kzg_instance.g_1.mul_bigint(Fr::from(42).into_bigint());

        assert_eq!(commit_result, expected_commit_result);
    }

    #[test]
    fn test_open() {
        let poly_evals = &[
            Fr::from(0),
            Fr::from(4),
            Fr::from(0),
            Fr::from(4),
            Fr::from(0),
            Fr::from(4),
            Fr::from(3),
            Fr::from(7),
        ];
        let unenc_taus = vec![Fr::from(5), Fr::from(2), Fr::from(3)];
        let kzg_instance = KZG::new(&MultilinearPoly::new(poly_evals.to_vec()), unenc_taus);

        let opening_values = &[Fr::from(6), Fr::from(4), Fr::from(0)];

        let open_result = kzg_instance.open(opening_values);

        let expected_open_result = Fr::from(72);

        assert_eq!(open_result, expected_open_result);
    }

    #[test]
    fn test_get_proof() {
        let poly_evals = &[
            Fr::from(0),
            Fr::from(4),
            Fr::from(0),
            Fr::from(4),
            Fr::from(0),
            Fr::from(4),
            Fr::from(3),
            Fr::from(7),
        ];
        let unenc_taus = vec![Fr::from(5), Fr::from(2), Fr::from(3)];
        let kzg_instance = KZG::new(&MultilinearPoly::new(poly_evals.to_vec()), unenc_taus);

        let opening_values = &[Fr::from(6), Fr::from(4), Fr::from(0)];

        let opened_value = kzg_instance.open(opening_values);

        let proof_result = kzg_instance.get_proof(opened_value, opening_values);

        let expected_quotients = vec![Fr::from(6), Fr::from(18), Fr::from(4)];

        let expected_quotients_g1 = kzg_instance
            .g_1
            .batch_mul(&expected_quotients)
            .into_iter()
            .map(|point| point.into_group())
            .collect();

        let expected_proof = Proof {
            quotients: expected_quotients_g1,
        };

        assert_eq!(proof_result.quotients, expected_proof.quotients);
    }

    #[test]
    fn test_verify() {
        let poly_evals = &[
            Fr::from(0),
            Fr::from(4),
            Fr::from(0),
            Fr::from(4),
            Fr::from(0),
            Fr::from(4),
            Fr::from(3),
            Fr::from(7),
        ];

        let unenc_taus = vec![Fr::from(5), Fr::from(2), Fr::from(3)];
        let kzg_instance = KZG::new(&MultilinearPoly::new(poly_evals.to_vec()), unenc_taus);
        let opening_values = &[Fr::from(6), Fr::from(4), Fr::from(0)];
        let commitment = kzg_instance.commit();
        let opened_value = kzg_instance.open(opening_values);
        let proof = kzg_instance.get_proof(opened_value, opening_values);

        let is_verified = kzg_instance.verify(commitment, opened_value, proof, opening_values);

        assert_eq!(is_verified, true);
    }

    #[test]
    fn test_dont_verify_invalid_proof() {
        let poly_evals = &[
            Fr::from(0),
            Fr::from(4),
            Fr::from(0),
            Fr::from(4),
            Fr::from(0),
            Fr::from(4),
            Fr::from(3),
            Fr::from(7),
        ];
        let unenc_taus = vec![Fr::from(5), Fr::from(2), Fr::from(3)];
        let kzg_instance = KZG::new(&MultilinearPoly::new(poly_evals.to_vec()), unenc_taus);
        let opening_values = &[Fr::from(6), Fr::from(4), Fr::from(0)];
        let commitment = kzg_instance.commit();
        let opened_value = kzg_instance.open(opening_values);

        let invalid_proof = Proof {
            quotients: vec![G1::generator(), G1::generator(), G1::generator()],
        };

        let is_verified =
            kzg_instance.verify(commitment, opened_value, invalid_proof, opening_values);

        assert_eq!(is_verified, false);
    }
}
