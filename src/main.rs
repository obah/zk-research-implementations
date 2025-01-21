use ark_bn254::Fq;
use implementations::univariate_polynomial::UnivariatePoly;

mod implementations;

fn check_test(x: Fq, poly: &UnivariatePoly<Fq>) {
    assert_eq!(
        poly.evaluate(x),
        (poly.evaluate(x - Fq::from(1)) + poly.evaluate(x - Fq::from(2)))
    );
}

fn fibonacci_check() {
    let fib_points = vec![
        Fq::from(1),
        Fq::from(1),
        Fq::from(2),
        Fq::from(3),
        Fq::from(5),
        Fq::from(8),
        Fq::from(13),
        Fq::from(21),
    ];

    let poly_points = fib_points
        .iter()
        .enumerate()
        .map(|(i, point)| (Fq::from(i as u64), *point))
        .collect();

    let polynomial = UnivariatePoly::interpolate(poly_points);

    check_test(Fq::from(2), &polynomial);
    check_test(Fq::from(5), &polynomial);
    check_test(Fq::from(7), &polynomial);
}

fn main() {
    println!("ZK Bootcamp");

    fibonacci_check();
}
