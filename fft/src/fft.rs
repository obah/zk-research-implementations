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
