use ark_ff::PrimeField;

enum Operation {
    Add,
    Mul,
}

struct Circuit {
    layers: Vec<Vec<Operation>>,
}

impl Circuit {
    fn new(layers: Vec<Vec<Operation>>) -> Self {
        Self { layers }
    }

    fn evaluate<F: PrimeField>(&self, inputs: Vec<F>) -> Vec<Vec<F>> {
        let mut current_inputs = inputs;
        let mut results = Vec::with_capacity(self.layers.len());

        for layer in &self.layers {
            assert!(
                current_inputs.len() >= 2 * layer.len(),
                "Insufficient inputs for layer with {} operations",
                layer.len()
            );

            let outputs: Vec<F> = current_inputs
                .chunks_exact(2)
                .zip(layer)
                .map(|(chunk, op)| match op {
                    Operation::Add => chunk[0] + chunk[1],
                    Operation::Mul => chunk[0] * chunk[1],
                })
                .collect();

            results.push(outputs.clone());
            current_inputs = outputs;
        }

        results
    }
}

#[cfg(test)]
mod test {
    use super::{Circuit, Operation};
    use ark_bn254::Fq;

    #[test]
    fn it_evaluates_the_circuit_correctly() {
        let structure: Vec<Vec<Operation>> = vec![
            vec![
                Operation::Mul,
                Operation::Mul,
                Operation::Mul,
                Operation::Mul,
            ],
            vec![Operation::Add, Operation::Add],
            vec![Operation::Add],
        ];

        let inputs: Vec<Fq> = vec![
            Fq::from(5),
            Fq::from(2),
            Fq::from(2),
            Fq::from(4),
            Fq::from(10),
            Fq::from(0),
            Fq::from(3),
            Fq::from(3),
        ];

        let expected_evaluations: Vec<Vec<Fq>> = vec![
            vec![Fq::from(10), Fq::from(8), Fq::from(0), Fq::from(9)],
            vec![Fq::from(18), Fq::from(9)],
            vec![Fq::from(27)],
        ];

        let circuit = Circuit::new(structure);

        let evaluations = circuit.evaluate(inputs);

        assert_eq!(evaluations, expected_evaluations);
    }
}
