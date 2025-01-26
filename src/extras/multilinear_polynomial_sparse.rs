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
