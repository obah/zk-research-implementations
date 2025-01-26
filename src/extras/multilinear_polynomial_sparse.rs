//implement sparse method
//implement desnse method
//implement partial and complete evaluation for both

use ark_ff::PrimeField;

//TODO use generic T to represent tuple of u8 in varying sizes instead of size 3 tuple - (F,F,F)
type Coefficients<F> = Vec<((u8, u8, u8), F)>;

struct MultilinearPolySparse<F: PrimeField> {
    coefficient: Coefficients<F>,
}

trait MultilinearPolySparseTrait<F> {
    fn new(coeff: Coefficients<F>) -> Self;

    fn evaluate(&self, points: (F, F, F)) -> F;

    fn partial_evaluate(&self, points: (F, F, F)) -> F;

    fn scalar_mul(&self, scalar: F) -> Self;

    fn interpolate(points: Coefficients<F>) -> Self;
}

impl<F: PrimeField> MultilinearPolySparseTrait<F> for MultilinearPolySparse<F> {
    fn new(coeff: Coefficients<F>) -> Self {
        MultilinearPolySparse { coefficient: coeff }
    }

    fn evaluate(&self, points: (F, F, F)) -> F {
        let (x, y, z) = points;

        self.coefficient
            .iter()
            .map(|(variables, coeff)| {
                let (a, b, c) = variables;

                *coeff * (x.pow([*a as u64])) * (y.pow([*b as u64])) * (z.pow([*c as u64]))
            })
            .sum()
    }

    fn partial_evaluate(&self, points: (F, F, F)) -> F {
        todo!()
    }

    fn scalar_mul(&self, scalar: F) -> Self {
        todo!()
    }

    fn interpolate(points: Coefficients<F>) -> Self {
        todo!()
    }
}

//? boolean hypercube from multilinear evaluation
// fn generate_bhc(bits: usize, poly_evaluation: Vec<F>) -> Self {
//     let size = 1 << bits;
//     // let mut binary_values = Vec::with_capacity(size);

//     // for i in 0..size {
//     //     let mut point = Vec::with_capacity(bits);

//     //     for j in (0..bits).rev() {
//     //         point.push(((i >> j) & 1) as u8);
//     //     }

//     //     binary_values.push(point);
//     // }

//     let binary_values: Vec<Vec<u8>> = (0..size)
//         .map(|i| (0..bits).rev().map(|j| ((i >> j) & 1) as u8).collect())
//         .collect();

//     // let hypercube: Hypercube<F> = binary_values
//     //     .iter()
//     //     .enumerate()
//     //     .map(|(i, point)| (point.clone(), poly_evaluation[i].clone()))
//     //     .collect();

//     let hypercube: Hypercube<F> = binary_values
//         .iter()
//         .enumerate()
//         .map(|(i, point)| (point.clone(), poly_evaluation[i])) // Avoid cloning F
//         .collect();

//     Self::new(hypercube)
// }
