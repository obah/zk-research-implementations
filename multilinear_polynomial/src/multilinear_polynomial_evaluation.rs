use ark_ff::PrimeField;
use std::ops::{Add, Mul, Sub};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Operation {
    Add,
    Mul,
}

impl Operation {
    pub fn apply<F: PrimeField>(self, a: F, b: F) -> F {
        match self {
            Operation::Add => a + b,
            Operation::Mul => a * b,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct MultilinearPoly<F: PrimeField> {
    pub evaluation: Vec<F>,
    pub num_of_vars: usize,
}

impl<F: PrimeField> MultilinearPoly<F> {
    pub fn new(evaluations: Vec<F>) -> Self {
        let num_of_vars: usize = evaluations.len().ilog2() as usize;

        if evaluations.len() != 1 << num_of_vars {
            panic!("Invalid evaluations");
        }

        Self {
            evaluation: evaluations,
            num_of_vars,
        }
    }

    fn pair_points(bit: usize, num_of_vars: usize) -> Vec<(usize, usize)> {
        let mut result = vec![];
        let target_hc = num_of_vars - 1;

        for val in 0..(1 << target_hc) {
            let inverted_index = num_of_vars - bit - 1;
            let insert_zero = insert_bit(val, inverted_index);
            let insert_one = insert_zero | (1 << inverted_index);
            result.push((insert_zero, insert_one));
        }
        result
    }

    pub fn partial_evaluate(&self, bit: usize, value: &F) -> Self {
        let mut result: Vec<F> = Vec::new();

        for (a, b) in MultilinearPoly::<F>::pair_points(bit, self.num_of_vars).into_iter() {
            let a = self.evaluation[a];
            let b = self.evaluation[b];

            result.push(a + *value * (b - a));
        }

        Self::new(result)
    }

    pub fn multi_partial_evaluate(&self, values: &[F]) -> Self {
        if values.len() > self.num_of_vars {
            panic!("Invalid number of values");
        }

        let mut poly = self.clone();

        for value in values {
            poly = poly.partial_evaluate(0, value);
        }

        poly
    }

    pub fn evaluate(&self, values: Vec<F>) -> F {
        if values.len() != self.num_of_vars {
            panic!("Invalid number of values");
        }

        let mut result = self.clone();

        for value in values.iter() {
            result = result.partial_evaluate(0, value);
        }

        result.evaluation[0]
    }

    pub fn scale(&self, value: F) -> Self {
        let result = self.evaluation.iter().map(|eval| *eval * value).collect();

        Self::new(result)
    }

    pub fn tensor_add_mul_polynomials(
        poly_a: &[F],
        poly_b: &[F],
        op: Operation,
    ) -> MultilinearPoly<F> {
        let new_eval: Vec<F> = poly_a
            .iter()
            .flat_map(|a| poly_b.iter().map(move |b| op.apply(*a, *b)))
            .collect();

        MultilinearPoly::new(new_eval)
    }
}

impl<F: PrimeField> Add for MultilinearPoly<F> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let result = self
            .evaluation
            .iter()
            .zip(other.evaluation.iter())
            .map(|(a, b)| *a + *b)
            .collect();

        MultilinearPoly::new(result)
    }
}

impl<F: PrimeField> Mul for MultilinearPoly<F> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let result = self
            .evaluation
            .iter()
            .zip(other.evaluation.iter())
            .map(|(a, b)| *a * *b)
            .collect();

        MultilinearPoly::new(result)
    }
}

impl<F: PrimeField> Sub for MultilinearPoly<F> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let result = self
            .evaluation
            .iter()
            .zip(other.evaluation.iter())
            .map(|(a, b)| *a - *b)
            .collect();

        MultilinearPoly::new(result)
    }
}

fn insert_bit(value: usize, bit: usize) -> usize {
    let high = value >> bit;
    let mask = (1 << bit) - 1;
    let low = value & mask;

    high << (bit + 1) | low
}

#[cfg(test)]
mod test {
    use super::*;
    use field_tracker::{print_summary, Ft};

    type Fq = Ft!(ark_bn254::Fq);

    #[test]
    fn it_partially_evaluates_any_multilinear() {
        let evaluations = vec![Fq::from(0), Fq::from(0), Fq::from(3), Fq::from(10)];
        let polynomial = MultilinearPoly::new(evaluations);

        let value_a = Fq::from(5);
        let bit_a = 0;

        let result = polynomial.partial_evaluate(bit_a, &value_a);

        assert_eq!(result.evaluation, vec![Fq::from(15), Fq::from(50)]);

        print_summary!();
    }

    #[test]
    fn it_fully_evaluates_any_multilinear() {
        let evaluations = vec![Fq::from(0), Fq::from(0), Fq::from(3), Fq::from(10)];
        let polynomial = MultilinearPoly::new(evaluations);

        let values = vec![Fq::from(5), Fq::from(1)];

        let result = polynomial.evaluate(values);

        assert_eq!(result, Fq::from(50));
    }
}
