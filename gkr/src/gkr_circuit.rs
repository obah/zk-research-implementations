use ark_ff::PrimeField;

#[derive(Clone, PartialEq)]
enum Operation {
    Add,
    Mul,
}

struct Gate<F: PrimeField> {
    l_input: F,
    r_input: F,
    output: F,
    ops: Operation,
}

struct Layer<F: PrimeField> {
    gates: Vec<Gate<F>>,
}

struct Circuit<F: PrimeField> {
    layers: Vec<Layer<F>>,
}

impl<F: PrimeField> Gate<F> {
    fn new(l_input: F, r_input: F, ops: Operation) -> Self {
        let output: F;

        match ops {
            Operation::Add => output = r_input + l_input,
            Operation::Mul => output = r_input * l_input,
        }

        Self {
            l_input,
            r_input,
            output,
            ops,
        }
    }
}

impl<F: PrimeField> Layer<F> {
    fn new(gates: Vec<Gate<F>>) -> Self {
        Self { gates }
    }

    fn evaluate(&self) -> Vec<F> {
        let result = self
            .gates
            .iter()
            .map(|gate| gate.output)
            .collect::<Vec<_>>();

        result
    }

    fn get_w_poly(&self) -> Vec<F> {
        self.evaluate()
    }

    fn get_add_mul_i(&self, op: Operation) -> Vec<u8> {
        let n_bits = self.get_bits_for_gates();
        let layer_size = 1 << n_bits;

        let mut poly_eval: Vec<u8> = vec![0; layer_size];

        let gate_values = self.gate_to_bits();

        for (i, gate) in self.gates.iter().enumerate() {
            if gate.ops == op {
                let gate_value = gate_values[i];
                poly_eval[gate_value] = 1;
            }
        }

        poly_eval
    }

    fn get_bits_for_gates(&self) -> u32 {
        let n_gates = self.gates.len();

        assert!(n_gates > 0, "There must be at least one gate in the layer.");

        let n_gates_log = n_gates.ilog2();
        let n_bits = n_gates_log + 1;

        match n_gates {
            1 => 3, // didnt recalucalate this since I know its always this.
            _ => n_gates_log + (n_bits * 2),
        }
    }

    fn gate_to_bits(&self) -> Vec<usize> {
        self.gates
            .iter()
            .enumerate()
            .map(|(i, _)| {
                let mut gate_value = 0;
                // output
                gate_value += i;

                let gate_id = i * 2;

                //inputs
                gate_value += 2 * gate_id + 1;

                gate_value
            })
            .collect()
    }
}

impl<F: PrimeField> Circuit<F> {
    fn new(structure: Vec<Vec<Operation>>) -> Self {
        let new_circuit = structure
            .iter()
            .map(|layer| {
                Layer::new(
                    layer
                        .iter()
                        .map(|gate| Gate::new(F::zero(), F::zero(), gate.clone()))
                        .collect(),
                )
            })
            .collect::<Vec<_>>();

        Self {
            layers: new_circuit,
        }
    }

    fn evaluate(&mut self, inputs: Vec<F>) -> Vec<Vec<F>> {
        let mut result = Vec::new();
        let mut current_inputs = inputs;

        for layer in self.layers.iter_mut() {
            for (i, gate) in layer.gates.iter_mut().enumerate() {
                let input_id = i * 2;

                gate.l_input = current_inputs[input_id];
                gate.r_input = current_inputs[input_id + 1];

                gate.output = match gate.ops {
                    Operation::Add => gate.l_input + gate.r_input,
                    Operation::Mul => gate.l_input * gate.r_input,
                };
            }

            let layer_outputs = layer.evaluate();
            result.push(layer_outputs.clone());
            current_inputs = layer_outputs;
        }
        result
    }
}

// use ark_ff::PrimeField;

// enum Operation {
//     Add,
//     Mul,
// }

// struct Circuit {
//     layers: Vec<Vec<Operation>>,
// }

// impl Circuit {
//     fn new(layers: Vec<Vec<Operation>>) -> Self {
//         Self { layers }
//     }

//     fn evaluate<F: PrimeField>(&self, inputs: Vec<F>) -> Vec<Vec<F>> {
//         let mut current_inputs = inputs;
//         let mut results = Vec::with_capacity(self.layers.len());

//         for layer in &self.layers {
//             assert!(
//                 current_inputs.len() >= 2 * layer.len(),
//                 "Insufficient inputs for layer with {} operations",
//                 layer.len()
//             );

//             let outputs: Vec<F> = current_inputs
//                 .chunks_exact(2)
//                 .zip(layer)
//                 .map(|(chunk, op)| match op {
//                     Operation::Add => chunk[0] + chunk[1],
//                     Operation::Mul => chunk[0] * chunk[1],
//                 })
//                 .collect();

//             results.push(outputs.clone());
//             current_inputs = outputs;
//         }

//         results
//     }
// }

#[cfg(test)]
mod test {
    use super::{Circuit, Gate, Layer, Operation};
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

        let mut circuit = Circuit::new(structure);

        let evaluations = circuit.evaluate(inputs);

        assert_eq!(evaluations, expected_evaluations);
    }

    #[test]
    fn it_returns_right_w_polys_for_each_layer() {
        let gate_1 = Gate::new(Fq::from(1), Fq::from(2), Operation::Add);
        let gate_2 = Gate::new(Fq::from(3), Fq::from(4), Operation::Mul);
        let gate_3 = Gate::new(Fq::from(5), Fq::from(6), Operation::Add);
        let gate_4 = Gate::new(Fq::from(7), Fq::from(8), Operation::Mul);

        let layer = Layer::new(vec![gate_1, gate_2, gate_3, gate_4]);

        let expected_w_poly = vec![Fq::from(3), Fq::from(12), Fq::from(11), Fq::from(56)];

        let w_poly = layer.get_w_poly();

        assert_eq!(w_poly, expected_w_poly);
    }

    // #[test]
    // fn it_returns_right_add_i_polys_for_each_layer() {
    //     let gate_1 = Gate::new(Fq::from(1), Fq::from(2), Operation::Add);
    //     let gate_2 = Gate::new(Fq::from(3), Fq::from(4), Operation::Mul);
    //     let gate_3 = Gate::new(Fq::from(5), Fq::from(6), Operation::Add);
    //     let gate_4 = Gate::new(Fq::from(7), Fq::from(8), Operation::Mul);

    //     let layer = Layer::new(vec![gate_1, gate_2, gate_3, gate_4]);

    //     // let expected_add_i_poly = vec![];

    //     let add_i_poly = layer.get_add_mul_i(Operation::Add);

    //     assert_eq!(add_i_poly.iter().sum(), 2);
    // }

    // #[test]
    // fn it_returns_right_mul_i_polys_for_each_layer() {
    //     let gate_1 = Gate::new(Fq::from(1), Fq::from(2), Operation::Add);
    //     let gate_2 = Gate::new(Fq::from(3), Fq::from(4), Operation::Mul);
    //     let gate_3 = Gate::new(Fq::from(5), Fq::from(6), Operation::Add);
    //     let gate_4 = Gate::new(Fq::from(7), Fq::from(8), Operation::Mul);

    //     let layer = Layer::new(vec![gate_1, gate_2, gate_3, gate_4]);

    //     let mul_i_poly = layer.get_add_mul_i(Operation::Mul);

    //     assert_eq!(mul_i_poly.iter().sum(), 2);
    // }

    #[test]
    fn it_returns_right_add_i_polys_for_each_layer() {
        let gate_1 = Gate::new(Fq::from(1), Fq::from(2), Operation::Add);
        let gate_2 = Gate::new(Fq::from(1), Fq::from(2), Operation::Mul);

        let layer_1 = Layer::new(vec![gate_1]);
        let layer_2 = Layer::new(vec![gate_2]);

        let expected_add_1_poly = vec![0, 1, 0, 0, 0, 0, 0, 0];
        let expected_add_2_poly = vec![0, 0, 0, 0, 0, 0, 0, 0];

        let add_1_poly = layer_1.get_add_mul_i(Operation::Add);
        let add_2_poly = layer_2.get_add_mul_i(Operation::Add);

        assert_eq!(expected_add_1_poly, add_1_poly);
        assert_eq!(expected_add_2_poly, add_2_poly);
    }

    #[test]
    fn it_returns_right_mul_i_polys_for_each_layer() {
        let gate_1 = Gate::new(Fq::from(1), Fq::from(2), Operation::Add);
        let gate_2 = Gate::new(Fq::from(1), Fq::from(2), Operation::Mul);

        let layer_1 = Layer::new(vec![gate_1]);
        let layer_2 = Layer::new(vec![gate_2]);

        let expected_mul_1_poly = vec![0, 0, 0, 0, 0, 0, 0, 0];
        let expected_mul_2_poly = vec![0, 1, 0, 0, 0, 0, 0, 0];

        let mul_1_poly = layer_1.get_add_mul_i(Operation::Mul);
        let mul_2_poly = layer_2.get_add_mul_i(Operation::Mul);

        assert_eq!(expected_mul_1_poly, mul_1_poly);
        assert_eq!(expected_mul_2_poly, mul_2_poly);
    }
}
