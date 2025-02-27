use ark_bls12_381::{Bls12_381, Fr, G1Projective as G1, G2Projective as G2};
use ark_ec::{pairing::Pairing, AffineRepr, PrimeGroup, ScalarMul};
use ark_ff::{PrimeField, UniformRand};
use gkr::{gkr_circuit::Operation, gkr_protocol::tensor_add_mul_polynomials};
use multilinear_polynomial::multilinear_polynomial_evaluation::MultilinearPoly;

struct Proof {
    quotients: Vec<G1>,
}

#[derive(Debug)]
struct LagrangeBasis {
    values: Vec<G1>,
    num_of_vars: usize,
}

struct KZG {
    polynomial: MultilinearPoly<Fr>,
    g_1: G1,
    g_2: G2,
    g1_taus: Vec<G1>,
    g2_taus: Vec<G2>,
    lagrange_basis: Vec<LagrangeBasis>,
}

impl KZG {
    fn new(polynomial: &MultilinearPoly<Fr>) -> Self {
        let g_1 = G1::generator();
        let g_2 = G2::generator();

        let (g1_taus, g2_taus, lagrange_basis) = KZG::run_trusted_setup(polynomial, g_1, g_2);

        Self {
            polynomial: polynomial.clone(),
            g_1,
            g_2,
            g1_taus,
            g2_taus,
            lagrange_basis,
        }
    }

    fn run_trusted_setup(
        poly: &MultilinearPoly<Fr>,
        g_1: G1,
        g_2: G2,
    ) -> (Vec<G1>, Vec<G2>, Vec<LagrangeBasis>) {
        let mut random_nums = Vec::with_capacity(poly.num_of_vars);

        let mut rng = ark_std::test_rng();

        for _ in 0..poly.num_of_vars {
            let random_scalar = Fr::rand(&mut rng);
            random_nums.push(random_scalar);
        }

        let g1_taus_affine = g_1.batch_mul(&random_nums);
        let g2_taus_affine = g_2.batch_mul(&random_nums);

        let g1_taus: Vec<G1> = g1_taus_affine
            .into_iter()
            .map(|point| point.into_group())
            .collect();

        let g2_taus: Vec<G2> = g2_taus_affine
            .into_iter()
            .map(|point| point.into_group())
            .collect();

        let mut all_lagrange_basis = Vec::new();

        let mut n_vars = poly.num_of_vars;
        while n_vars > 0 {
            let evaluation_g1 = get_lagrange_basis(n_vars, &random_nums, g_1);

            let lagrange_basis = LagrangeBasis {
                values: evaluation_g1,
                num_of_vars: n_vars,
            };

            all_lagrange_basis.push(lagrange_basis);

            n_vars -= 1;
        }

        (g1_taus, g2_taus, all_lagrange_basis)
    }

    fn commit(&self) -> G1 {
        evaluate_poly_with_l_basis_in_g1(&self.polynomial.evaluation, &self.lagrange_basis[0])
    }

    fn open(&self, opening_values: &[Fr]) -> Fr {
        // let result = self.polynomial.evaluate(opening_values.to_vec());
        self.polynomial.evaluate(opening_values.to_vec())

        // let result_g1 = self.g_1.mul_bigint(result.into_bigint());

        // result_g1
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

        let mut num_of_vars = opening_values.len() - 1;
        for value in opening_values {
            let quotient = get_quotient(&poly_minus_v, 0);

            let quotient_eval;

            if num_of_vars > 0 {
                let lagrange_basis = self
                    .lagrange_basis
                    .iter()
                    .find(|basis| basis.num_of_vars == num_of_vars)
                    .unwrap();

                quotient_eval = evaluate_poly_with_l_basis_in_g1(&quotient, lagrange_basis);

                num_of_vars -= 1;
            } else {
                quotient_eval = self.g_1.mul_bigint(quotient[0].into_bigint());
                // quotient_eval = encrypt_value(quotient[0], self.g_1);
            }

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

        // let lhs = commitment - encrypt_value(opened_value, self.g_1);
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

fn evaluate_poly_with_l_basis_in_g1(poly_evaluations: &[Fr], lagrange_basis: &LagrangeBasis) -> G1 {
    if poly_evaluations.len() != lagrange_basis.values.len() {
        panic!("invalid polynomial or lagrange basis");
    }

    poly_evaluations
        .iter()
        .zip(lagrange_basis.values.iter())
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

    tensor_add_mul_polynomials(&blow_up_poly, poly, Operation::Mul).evaluation
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

// fn encrypt_value(value: Fr, g_1: G1) -> G1 {
//     g_1.mul_bigint(value.into_bigint())
// }

#[cfg(test)]
mod test {
    use ark_bls12_381::{Fr, G1Projective};
    use ark_ec::{PrimeGroup, ScalarMul};
    use ark_ff::PrimeField;

    use crate::kzg_pcs::kzg::LagrangeBasis;

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

        let lagrange_basis_evals = get_lagrange_basis(n_vars, unenc_taus, g_1);

        let lagrange_basis = LagrangeBasis {
            values: lagrange_basis_evals,
            num_of_vars: n_vars,
        };

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

        let kzg_instance = KZG::new(&MultilinearPoly::new(poly_evals.to_vec()));

        let commit_result = kzg_instance.commit();

        let expected_commit_result = kzg_instance.g_1.mul_bigint(Fr::from(42).into_bigint());

        assert_eq!(commit_result, expected_commit_result);
        ////! cant test yet becuase of the randomness
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

        let kzg_instance = KZG::new(&MultilinearPoly::new(poly_evals.to_vec()));

        let opening_values = &[Fr::from(6), Fr::from(4), Fr::from(0)];

        let open_result = kzg_instance.open(opening_values);

        // let expected_open_result = kzg_instance.g_1.mul_bigint(Fr::from(72).into_bigint());
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

        let kzg_instance = KZG::new(&MultilinearPoly::new(poly_evals.to_vec()));

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

        ////! cant test yet becuase of the randomness
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

        let kzg_instance = KZG::new(&MultilinearPoly::new(poly_evals.to_vec()));
        let opening_values = &[Fr::from(6), Fr::from(4), Fr::from(0)];
        let commitment = kzg_instance.commit();
        let opened_value = kzg_instance.open(opening_values);
        let proof = kzg_instance.get_proof(opened_value, opening_values);

        let is_verified = kzg_instance.verify(commitment, opened_value, proof, opening_values);

        assert_eq!(is_verified, true);

        ////! cant test yet becuase of the randomness
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

        let kzg_instance = KZG::new(&MultilinearPoly::new(poly_evals.to_vec()));
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
