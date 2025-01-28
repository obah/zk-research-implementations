use ark_ff::PrimeField;
use std::ops::{Add, AddAssign};

type Hypercube<F> = Vec<(Vec<u8>, F)>;

struct BHC<F> {
    bits: Hypercube<F>,
}

#[derive(Clone)]
pub struct MultilinearPoly<F: PrimeField> {
    pub evaluation: Vec<F>,
}

impl<F: PrimeField> BHC<F> {
    fn new(hypercube: Hypercube<F>) -> Self {
        BHC { bits: hypercube }
    }

    fn generate_bhc(poly_evaluation: Vec<F>) -> Self {
        let bits = poly_evaluation.len().ilog2() as usize;
        let size = 1 << bits;

        let hypercube: Hypercube<F> = (0..size)
            .map(|i| {
                let point = (0..bits).rev().map(|j| ((i >> j) & 1) as u8).collect();
                (point, poly_evaluation[i])
            })
            .collect();

        Self::new(hypercube)
    }

    fn pair_points(&self, bit: u8) -> Vec<(F, F)> {
        let mut pairs = Vec::new();
        let pair_index = 1 << bit;

        for i in 0..pair_index {
            if i + pair_index < self.bits.len() {
                let (_, a) = &self.bits[i];
                let (_, b) = &self.bits[i + pair_index];
                pairs.push((a.clone(), b.clone()));
            }
        }

        pairs
    }
}

impl<F: PrimeField> MultilinearPoly<F> {
    pub fn new(evaluations: Vec<F>) -> Self {
        MultilinearPoly {
            evaluation: evaluations,
        }
    }

    fn iterpolate(points: (F, F), value: F) -> F {
        let (y_0, y_1) = points;

        y_0 + (value * (y_1 - y_0))
    }

    pub fn partial_evaluate(&self, value: F, bit: u8) -> Self {
        let bhc = BHC::generate_bhc(self.evaluation.clone());

        let paired_evaluations = bhc
            .pair_points(bit)
            .iter()
            .map(|point| Self::iterpolate(*point, value))
            .collect();

        Self::new(paired_evaluations)
    }

    // pub fn evaluate(&self, values: Vec<F>) -> F {
    //     let mut result = self;

    //     let mut bits = result.evaluation.len().ilog2() - 1;

    //     for value in values.iter() {
    //         result = result.partial_evaluate(*value, bits.try_into().unwrap());

    //         if bits == 0 {
    //             break;
    //         } else {
    //             bits -= 1;
    //         }
    //     }

    //     result.evaluation[0]
    // }

    pub fn evaluate(&self, values: Vec<F>) -> F {
        let mut result = self.clone();

        let mut bits = result.evaluation.len().ilog2() - 1;

        for value in values.iter() {
            result = result.partial_evaluate(*value, bits.try_into().unwrap());

            if bits == 0 {
                break;
            } else {
                bits -= 1;
            }
        }

        // Return the final evaluation result
        result.evaluation[0]
    }
}

impl<F: PrimeField> Add for MultilinearPoly<F> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let mut result = vec![F::zero(); self.evaluation.len().max(other.evaluation.len())];

        for (i, &value) in self.evaluation.iter().enumerate() {
            result[i] += value;
        }

        for (i, &value) in other.evaluation.iter().enumerate() {
            result[i] += value;
        }

        MultilinearPoly::new(result)
    }
}

impl<F: PrimeField> AddAssign for MultilinearPoly<F> {
    fn add_assign(&mut self, other: Self) {
        let mut result = vec![F::zero(); self.evaluation.len().max(other.evaluation.len())];

        for (i, &value) in self.evaluation.iter().enumerate() {
            result[i] += value;
        }

        for (i, &value) in other.evaluation.iter().enumerate() {
            result[i] += value;
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bn254::Fq;

    #[test]
    fn it_pair_points_correctly() {
        let evaluations = vec![Fq::from(0), Fq::from(1), Fq::from(2), Fq::from(3)];
        let bhc = BHC::generate_bhc(evaluations);

        let pairs = bhc.pair_points(1);
        assert_eq!(
            pairs,
            vec![(Fq::from(0), Fq::from(2)), (Fq::from(1), Fq::from(3))]
        );
    }

    #[test]
    fn it_partially_evaluates_any_multilinear() {
        let evaluations = vec![Fq::from(0), Fq::from(0), Fq::from(3), Fq::from(10)];
        let polynomial = MultilinearPoly::new(evaluations);

        let value_a = Fq::from(5);
        let bit_a = 1;

        let result = polynomial.partial_evaluate(value_a, bit_a);

        assert_eq!(result.evaluation, vec![Fq::from(15), Fq::from(50)]);
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
