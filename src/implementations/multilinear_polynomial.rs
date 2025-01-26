use ark_ff::PrimeField;

type Hypercube<F> = Vec<(Vec<u8>, F)>;

struct BHC<F> {
    bits: Hypercube<F>,
}

struct MultilinearPoly<F: PrimeField> {
    evaluation: Vec<F>,
}

impl<F: PrimeField> BHC<F> {
    fn new(hypercube: Hypercube<F>) -> Self {
        BHC { bits: hypercube }
    }

    fn generate_bhc(bits: usize, poly_evaluation: Vec<F>) -> Self {
        let size = 1 << bits;
        let mut binary_values = Vec::with_capacity(size);

        for i in 0..size {
            let mut point = Vec::with_capacity(bits);

            for j in (0..bits).rev() {
                point.push(((i >> j) & 1) as u8);
            }

            binary_values.push(point);
        }

        let hypercube: Hypercube<F> = binary_values
            .iter()
            .enumerate()
            .map(|(i, point)| (point.clone(), poly_evaluation[i].clone()))
            .collect();

        Self::new(hypercube)
    }

    fn pair_points(&mut self, bit: u8) -> Vec<(F, F)> {
        let mut pairs = Vec::new();
        let mut pair_index = 1 << bit;

        while pair_index > 0 && self.bits.len() > 1 {
            if pair_index < self.bits.len() {
                let (_, a) = self.bits.remove(0);
                let (_, b) = self.bits.remove(pair_index - 1);
                pairs.push((a, b));
            } else {
                break;
            }
            pair_index -= 1;
        }

        pairs
    }
}

impl<F: PrimeField> MultilinearPoly<F> {
    fn new(evaluations: Vec<F>) -> Self {
        MultilinearPoly {
            evaluation: evaluations,
        }
    }

    fn iterpolate(points: (F, F), value: F) -> F {
        let (y_0, y_1) = points;

        y_0 + (value * (y_1 - y_0))
    }

    fn partial_evaluate(&self, value: F, bit: u8) -> Self {
        let bits = self.evaluation.len().ilog2() as usize;

        let mut bhc = BHC::generate_bhc(bits, self.evaluation.clone());

        let paired_evaluations = bhc
            .pair_points(bit)
            .iter()
            .map(|point| Self::iterpolate(*point, value))
            .collect();

        Self::new(paired_evaluations)
    }

    fn evaluate(self, values: Vec<F>) -> F {
        let mut result = self;

        let mut bits = result.evaluation.len().ilog2() - 1;

        for value in values.iter() {
            result = result.partial_evaluate(*value, bits.try_into().unwrap());

            if bits == 0 {
                break;
            } else {
                bits -= 1;
            }
        }

        result.evaluation[0]
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bn254::Fq;

    #[test]
    fn it_pair_points_correctly() {
        let evaluations = vec![Fq::from(0), Fq::from(1), Fq::from(2), Fq::from(3)];
        let mut bhc = BHC::generate_bhc(2, evaluations);

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
