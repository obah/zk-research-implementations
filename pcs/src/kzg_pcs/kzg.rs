use ark_bls12_381::{Bls12_381, Fr, G1Projective as G1, G2Projective as G2};
use ark_ec::{pairing::Pairing, AffineRepr, PrimeGroup, ScalarMul};
use ark_ff::{PrimeField, UniformRand};
use gkr::{gkr_circuit::Operation, gkr_protocol::tensor_add_mul_polynomials};
use multilinear_polynomial::multilinear_polynomial_evaluation::MultilinearPoly;

struct Proof {
    quotients: Vec<G1>,
}

struct KZG {
    polynomial: MultilinearPoly<Fr>,
    g_1: G1,
    g_2: G2,
    g1_taus: Vec<G1>,
    g2_taus: Vec<G2>,
    lagrange_basis: Vec<G1>,
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
            g1_taus: g1_taus,
            g2_taus: g2_taus,
            lagrange_basis,
        }
    }

    fn run_trusted_setup(
        poly: &MultilinearPoly<Fr>,
        g_1: G1,
        g_2: G2,
    ) -> (Vec<G1>, Vec<G2>, Vec<G1>) {
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

        let lagrange_basis = get_lagrange_basis(&poly, &random_nums, g_1);

        (g1_taus, g2_taus, lagrange_basis)
    }

    fn commit(&self) -> G1 {
        self.polynomial
            .evaluation
            .iter()
            .zip(self.lagrange_basis.iter())
            .map(|(a, b)| b.mul_bigint(a.into_bigint()))
            .sum()
    }

    fn open(&self, opening_values: &[Fr]) -> G1 {
        let result = self.polynomial.evaluate(opening_values.to_vec());

        let result_g1 = self.g_1.mul_bigint(result.into_bigint());

        result_g1
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

        for value in opening_values {
            let quotient = get_quotient(&poly_minus_v, 0);

            //todo before pushing quotient,
            //blow up again to add the missing values
            //then mul with all taus to get the constant values

            q_i.push(quotient);

            let remainder_poly = get_remainder(&poly_minus_v, *value, 0);

            poly_minus_v = MultilinearPoly::new(remainder_poly);
        }

        Proof { quotients: q_i }
    }

    fn verify(
        &self,
        commitment: G1,
        opened_value: G1,
        proof: Proof,
        opening_values: &[Fr],
    ) -> bool {
        if proof.quotients.len() != opening_values.len() {
            panic!("num of quotients in proof not equal to num of opening values");
        }

        let lhs = commitment - opened_value;
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

fn get_remainder(poly: &MultilinearPoly<Fr>, value: Fr, bit: usize) -> Vec<Fr> {
    poly.partial_evaluate(bit, &value).evaluation
}

fn get_quotient(poly: &MultilinearPoly<Fr>, bit: usize) -> Vec<Fr> {
    let mut eval_0 = poly.partial_evaluate(bit, &Fr::from(0));
    let mut eval_1 = poly.partial_evaluate(bit, &Fr::from(1));
    let blown_up = false;

    if eval_0.evaluation.len() > eval_1.evaluation.len() {
        eval_1 = MultilinearPoly::new(blow_up_poly(&eval_1.evaluation));
    } else if eval_0.evaluation.len() < eval_1.evaluation.len() {
        eval_0 = MultilinearPoly::new(blow_up_poly(&eval_0.evaluation));
    }

    let mut result = (eval_1 - eval_0).evaluation;

    if blown_up {
        result = result
            .chunks(result.len() / 2)
            .map(|chunk| {
                assert!(chunk[0] == chunk[1], "Chunk elements are not equal!");
                chunk[0]
            })
            .collect();
    }

    result
}

fn blow_up_poly(poly: &[Fr]) -> Vec<Fr> {
    let blow_up_factor = vec![Fr::from(1), Fr::from(1)];

    tensor_add_mul_polynomials(poly, &blow_up_factor, Operation::Mul).evaluation
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

fn get_lagrange_basis(poly: &MultilinearPoly<Fr>, unenc_taus: &[Fr], g_1: G1) -> Vec<G1> {
    let mut result = Vec::with_capacity(1 << poly.num_of_vars);
    let bhc = generate_bhc(poly.num_of_vars);

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
    use ark_bls12_381::Fr;

    use super::blow_up_poly;

    #[test]
    fn test_blow_up_poly() {
        let poly = vec![Fr::from(0), Fr::from(4)];

        let expected_poly = vec![Fr::from(0), Fr::from(4), Fr::from(0), Fr::from(4)];

        let result = blow_up_poly(&poly);

        assert_eq!(expected_poly, result);
    }
}
