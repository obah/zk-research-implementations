use ark_bn254::Fq;
use gkr::{
    gkr_circuit::{Circuit, Operation},
    gkr_protocol::{prove, verify},
};
use multilinear_polynomial::multilinear_polynomial_evaluation::MultilinearPoly;

fn main() {
    let eval_1 = vec![Fq::from(0), Fq::from(3), Fq::from(2), Fq::from(5)];
    let eval_2 = vec![Fq::from(0), Fq::from(3), Fq::from(2), Fq::from(5)];
    let eval_3 = vec![Fq::from(0), Fq::from(3), Fq::from(2), Fq::from(5)];
    let eval_4 = vec![Fq::from(0), Fq::from(3), Fq::from(2), Fq::from(5)];

    let mul_p_1 = MultilinearPoly::new(eval_1);
    let mul_p_2 = MultilinearPoly::new(eval_2);
    let mul_p_3 = MultilinearPoly::new(eval_3);
    let mul_p_4 = MultilinearPoly::new(eval_4);

    let circuit_structure: Vec<Vec<Operation>> =
        vec![vec![Operation::Mul, Operation::Mul], vec![Operation::Add]];

    let inputs: Vec<Fq> = vec![Fq::from(1), Fq::from(2), Fq::from(3), Fq::from(4)];

    // let circuit_structure: Vec<Vec<Operation>> = vec![
    //     vec![
    //         Operation::Mul,
    //         Operation::Mul,
    //         Operation::Mul,
    //         Operation::Mul,
    //     ],
    //     vec![Operation::Add, Operation::Add],
    //     vec![Operation::Add],
    // ];

    // let inputs: Vec<Fq> = vec![
    //     Fq::from(5),
    //     Fq::from(2),
    //     Fq::from(2),
    //     Fq::from(4),
    //     Fq::from(10),
    //     Fq::from(0),
    //     Fq::from(3),
    //     Fq::from(3),
    // ];

    let mut circuit = Circuit::new(circuit_structure);

    let proof = prove(&mut circuit, &inputs);

    let is_verified = verify(proof, circuit, &inputs);

    println!("is verified: {is_verified}");
}
