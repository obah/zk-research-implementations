use ark_bn254::Fq;
use gkr::{
    gkr_circuit::{Circuit, Operation},
    gkr_protocol::prove,
};

fn main() {
    let circuit_structure: Vec<Vec<Operation>> =
        vec![vec![Operation::Mul, Operation::Mul], vec![Operation::Add]];

    let inputs: Vec<Fq> = vec![Fq::from(1), Fq::from(2), Fq::from(3), Fq::from(4)];

    let mut circuit = Circuit::new(circuit_structure);

    let proof = prove(&mut circuit, inputs);

    // let is_verified = verify(proof, circuit);
}
