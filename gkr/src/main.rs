use ark_bn254::Fq;
use gkr::{
    gkr_circuit::{Circuit, Operation},
    gkr_protocol::{prove, verify},
};

fn main() {
    // let circuit_structure: Vec<Vec<Operation>> =
    //     vec![vec![Operation::Mul, Operation::Mul], vec![Operation::Add]];

    // let inputs: Vec<Fq> = vec![Fq::from(1), Fq::from(2), Fq::from(3), Fq::from(4)];

    let circuit_structure: Vec<Vec<Operation>> = vec![
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

    let mut circuit = Circuit::new(circuit_structure);

    let proof = prove(&mut circuit, &inputs);

    let is_verified = verify(proof, circuit, &inputs);

    println!("is verified: {is_verified}");
}
