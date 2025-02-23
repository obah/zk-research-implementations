use ark_ff::PrimeField;
use multilinear_polynomial::multilinear_polynomial_evaluation::MultilinearPoly;

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

#[derive(Debug, Clone)]
pub struct Gate<F: PrimeField> {
    pub l_input: F,
    pub r_input: F,
    pub output: F,
    pub op: Operation,
}

impl<F: PrimeField> Gate<F> {
    pub fn new(l_input: F, r_input: F, op: Operation) -> Self {
        let output = op.apply(l_input, r_input);

        Self {
            l_input,
            r_input,
            output,
            op,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Layer<F: PrimeField> {
    pub gates: Vec<Gate<F>>,
}

impl<F: PrimeField> Layer<F> {
    pub fn new(gates: Vec<Gate<F>>) -> Self {
        Self { gates }
    }

    pub fn get_layer_poly(&self) -> Vec<F> {
        self.gates.iter().map(|gate| gate.output).collect()
    }

    pub fn get_add_mul_i(&self, op: Operation) -> MultilinearPoly<F> {
        let n_bits = self.get_bits_for_gates();
        let layer_size = 1 << n_bits;
        let mut poly_eval = vec![F::zero(); layer_size];

        let gate_values = self.gate_to_bits();
        for (gate_value, gate) in gate_values.into_iter().zip(&self.gates) {
            if gate.op == op {
                poly_eval[gate_value] = F::one();
            }
        }

        MultilinearPoly::new(poly_eval)
    }

    fn get_bits_for_gates(&self) -> u32 {
        let n_gates = self.gates.len();
        assert!(n_gates > 0, "There must be at least one gate in the layer.");

        if n_gates == 1 {
            3
        } else {
            let n_gates_log = n_gates.ilog2();
            let n_bits = n_gates_log + 1;
            n_gates_log + (n_bits * 2)
        }
    }

    fn gate_to_bits(&self) -> Vec<usize> {
        let n_gates = self.gates.len();
        let n_gates_log = n_gates.ilog2();

        let mut gate_decimal_values = Vec::new();

        for (idx, _) in self.gates.iter().enumerate() {
            let mut gate_binary_values = Vec::new();

            gate_binary_values.push(idx);
            gate_binary_values.push(2 * idx);
            gate_binary_values.push(2 * idx + 1);

            let segments: Vec<(usize, u32)> = gate_binary_values
                .iter()
                .enumerate()
                .map(|(i, val)| {
                    if n_gates == 1 {
                        (*val, 1)
                    } else {
                        if i == 0 {
                            (*val, n_gates_log)
                        } else {
                            (*val, n_gates_log + 1)
                        }
                    }
                })
                .collect();

            let decimal_value = segments
                .iter()
                .fold(0, |acc, &(value, width)| (acc << width) | value);

            gate_decimal_values.push(decimal_value);
        }

        gate_decimal_values
    }
}

#[derive(Debug, Clone)]
pub struct Circuit<F: PrimeField> {
    pub layers: Vec<Layer<F>>,
}

impl<F: PrimeField> Circuit<F> {
    pub fn new(structure: Vec<Vec<Operation>>) -> Self {
        let layers = structure
            .into_iter()
            .map(|ops_layer| {
                let gates = ops_layer
                    .into_iter()
                    .map(|op| Gate::new(F::zero(), F::zero(), op))
                    .collect();
                Layer::new(gates)
            })
            .collect();
        Self { layers }
    }

    pub fn evaluate(&mut self, inputs: &[F]) -> Vec<Vec<F>> {
        let mut result = Vec::new();
        let mut current_inputs = inputs.to_vec();

        for layer in &mut self.layers {
            for (gate, input_pair) in layer.gates.iter_mut().zip(current_inputs.chunks_exact(2)) {
                let (l_input, r_input) = (input_pair[0], input_pair[1]);
                gate.l_input = l_input;
                gate.r_input = r_input;
                gate.output = gate.op.apply(l_input, r_input);
            }
            let layer_outputs = layer.get_layer_poly();
            result.push(layer_outputs.clone());
            current_inputs = layer_outputs;
        }
        result
    }
}

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

        let evaluations = circuit.evaluate(&inputs);

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

        let w_poly = layer.get_layer_poly();

        assert_eq!(w_poly, expected_w_poly);
    }

    #[test]
    fn it_returns_right_add_i_polys_for_each_layer() {
        let gate_1 = Gate::new(Fq::from(1), Fq::from(2), Operation::Add);
        let gate_2 = Gate::new(Fq::from(1), Fq::from(2), Operation::Mul);

        let layer_1 = Layer::new(vec![gate_1]);
        let layer_2 = Layer::new(vec![gate_2]);

        let expected_add_1_poly = vec![
            Fq::from(0),
            Fq::from(1),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
        ];
        let expected_add_2_poly = vec![Fq::from(0); 8];

        let add_1_poly = layer_1.get_add_mul_i(Operation::Add);
        let add_2_poly = layer_2.get_add_mul_i(Operation::Add);

        assert_eq!(expected_add_1_poly, add_1_poly.evaluation);
        assert_eq!(expected_add_2_poly, add_2_poly.evaluation);
    }

    #[test]
    fn it_returns_right_mul_i_polys_for_each_layer() {
        let gate_1 = Gate::new(Fq::from(1), Fq::from(2), Operation::Add);
        let gate_2 = Gate::new(Fq::from(1), Fq::from(2), Operation::Mul);

        let layer_1 = Layer::new(vec![gate_1]);
        let layer_2 = Layer::new(vec![gate_2]);

        let expected_mul_1_poly = vec![Fq::from(0); 8];
        let expected_mul_2_poly = vec![
            Fq::from(0),
            Fq::from(1),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
            Fq::from(0),
        ];

        let mul_1_poly = layer_1.get_add_mul_i(Operation::Mul);
        let mul_2_poly = layer_2.get_add_mul_i(Operation::Mul);

        assert_eq!(expected_mul_1_poly, mul_1_poly.evaluation);
        assert_eq!(expected_mul_2_poly, mul_2_poly.evaluation);
    }
}
