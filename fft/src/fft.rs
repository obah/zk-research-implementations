use ark_ff::PrimeField;
use univariate_polynomial::univariate_polynomial_dense::UnivariatePoly;

//takes poly in coeff form
fn fft_evaluate<F: PrimeField>(poly: &UnivariatePoly<F>) -> Vec<F> {
    let n = poly.coefficient.len();

    if !n.is_power_of_two() {
        panic!("Length of coefficients must be a power of 2");
    }

    if n == 1 {
        return poly.coefficient.clone();
    }

    let roots_of_unity = F::get_root_of_unity(n.try_into().unwrap()).unwrap();

    let (p_even, p_odd) = split_poly(&poly.coefficient);

    let p_even_poly = UnivariatePoly::new(p_even);
    let p_odd_poly = UnivariatePoly::new(p_odd);

    let (y_even, y_odd) = (fft_evaluate(&p_even_poly), fft_evaluate(&p_odd_poly));

    let mut y_points = Vec::with_capacity(n);

    for j in 0..n / 2 {
        y_points[j] = y_even[j] + (roots_of_unity[j] * y_odd[j]);
        y_points[j + (n / 2)] = y_even[j] - (roots_of_unity[j] * y_odd[j]);
    }

    y_points
}

//takes poly in evaluation form
fn fft_interpolate<F: PrimeField>(evaluations: &[F]) -> UnivariatePoly<F> {
    let n = evaluations.len();

    if !n.is_power_of_two() {
        panic!("Length of coefficients must be a power of 2");
    }

    if n == 1 {
        return UnivariatePoly::new(evaluations.to_vec());
    }

    //todo change this based on algorithm
    let roots_of_unity = F::get_root_of_unity(n.try_into().unwrap()).unwrap();

    let (p_even, p_odd) = split_poly(evaluations);

    let (y_even, y_odd) = (fft_interpolate(&p_even), fft_interpolate(&p_odd));

    let mut y_points = Vec::with_capacity(n);

    for j in 0..n / 2 {
        y_points[j] = y_even[j] + (roots_of_unity[j] * y_odd[j]);
        y_points[j + (n / 2)] = y_even[j] - (roots_of_unity[j] * y_odd[j]);
    }

    UnivariatePoly::new(y_points)
}

fn split_poly<F: PrimeField>(poly: &[F]) -> (Vec<F>, Vec<F>) {
    let poly_even = poly
        .iter()
        .enumerate()
        .filter(|(index, _)| index % 2 == 0)
        .map(|(_, coeff)| *coeff)
        .collect();

    let poly_odd = poly
        .iter()
        .enumerate()
        .filter(|(index, _)| index % 2 == 1)
        .map(|(_, coeff)| *coeff)
        .collect();

    ((poly_even), (poly_odd))
}

#[cfg(test)]

mod test {
    use ark_bn254::Fq;
    use univariate_polynomial::univariate_polynomial_dense::UnivariatePoly;

    use super::{fft_evaluate, fft_interpolate, split_poly};

    #[test]
    fn it_splits_poly_correctly() {
        //test with P = x3 + 2x2 -14x + 2
        let poly = &[Fq::from(2), Fq::from(-14), Fq::from(2), Fq::from(1)];

        let (p_even, p_odd) = split_poly(&poly);

        let expected_p_even = vec![Fq::from(2), Fq::from(2)];

        let expected_p_odd = vec![Fq::from(-14), Fq::from(1)];

        assert_eq!(p_even, expected_p_even);
        assert_eq!(p_odd, expected_p_odd);
    }

    #[test]
    fn it_evaluates_poly() {
        //test with P = 4x3 + 3x2 + 2x + 1
        let points = vec![Fq::from(1), Fq::from(2), Fq::from(3), Fq::from(4)];
        let poly = UnivariatePoly::new(points);

        let evaluations = fft_evaluate(&poly);

        let expected_evaluations = vec![Fq::from(), Fq::from(), Fq::from(), Fq::from()];

        assert_eq!(evaluations, expected_evaluations);
    }

    #[test]
    fn it_interpolates_poly() {
        //same poly as evaluate test
        let evaluations = vec![Fq::from(), Fq::from(), Fq::from(), Fq::from()];

        let poly = fft_interpolate(&evaluations);

        let expected_points = vec![Fq::from(1), Fq::from(2), Fq::from(3), Fq::from(4)];

        assert_eq!(poly.coefficient, expected_points);
    }
}
