use ark_bn254::Fr;
use ark_ff::{FftField, Field};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use univariate_polynomial::univariate_polynomial_dense::UnivariatePoly;

fn dft(values: &[Fr], root: Fr) -> Vec<Fr> {
    let n = values.len();

    if n == 1 {
        return values.to_vec();
    }

    let (even, odd) = split_poly(values);

    let root_sq = root.pow([2]);

    let y_even = dft(&even, root_sq);
    let y_odd = dft(&odd, root_sq);

    let mut y = vec![Fr::from(0); n];

    for j in 0..n / 2 {
        let twiddle = root.pow([j as u64]);
        y[j] = y_even[j] + twiddle * y_odd[j];
        y[j + n / 2] = y_even[j] - twiddle * y_odd[j];
    }

    y
}

fn fft_evaluate(poly: &UnivariatePoly<Fr>) -> Vec<Fr> {
    let n = poly.coefficient.len();

    if !n.is_power_of_two() {
        panic!("Length must be a power of 2");
    }

    let omega = Fr::get_root_of_unity(n as u64).unwrap();

    dft(&poly.coefficient, omega)
}

fn fft_interpolate(evaluations: &[Fr]) -> UnivariatePoly<Fr> {
    let n = evaluations.len();

    if !n.is_power_of_two() {
        panic!("Length must be a power of 2");
    }

    let omega_inv = Fr::get_root_of_unity(n as u64).unwrap().inverse().unwrap();

    let inv_n = Fr::from(n as u64).inverse().unwrap();

    let mut coeffs = dft(evaluations, omega_inv);

    for coeff in &mut coeffs {
        *coeff *= inv_n;
    }
    UnivariatePoly::new(coeffs)
}

fn split_poly(poly: &[Fr]) -> (Vec<Fr>, Vec<Fr>) {
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

fn get_interpolation_roots(n: usize) -> Vec<Fr> {
    let domain = GeneralEvaluationDomain::<Fr>::new(n).unwrap();

    let omega_inv = domain.group_gen_inv();

    let inv_n = domain.size_as_field_element().inverse().unwrap();

    (0..n).map(|k| omega_inv.pow([k as u64]) * inv_n).collect()
}

#[cfg(test)]
mod test {
    use ark_bn254::Fr;
    use ark_ff::{FftField, Field};
    use univariate_polynomial::univariate_polynomial_dense::UnivariatePoly;

    use super::{fft_evaluate, fft_interpolate, split_poly};

    #[test]
    fn test_splits_poly_correctly() {
        //test with P = x3 + 2x2 -14x + 2
        let poly = &[Fr::from(2), Fr::from(-14), Fr::from(2), Fr::from(1)];

        let (p_even, p_odd) = split_poly(poly);

        let expected_p_even = vec![Fr::from(2), Fr::from(2)];

        let expected_p_odd = vec![Fr::from(-14), Fr::from(1)];

        assert_eq!(p_even, expected_p_even);
        assert_eq!(p_odd, expected_p_odd);
    }

    #[test]

    fn test_evaluates_poly() {
        let coefficients = vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4)];
        let poly = UnivariatePoly::new(coefficients);

        let evaluations = fft_evaluate(&poly);

        let n = 4;
        let omega = Fr::get_root_of_unity(n as u64).unwrap();
        let roots: Vec<Fr> = (0..n).map(|i| omega.pow(&[i as u64])).collect();

        let expected_evaluations: Vec<Fr> = roots
            .iter()
            .map(|x| {
                let x2 = x.square();
                let x3 = x2 * x;
                Fr::from(1) + Fr::from(2) * x + Fr::from(3) * x2 + Fr::from(4) * x3
            })
            .collect();

        assert_eq!(evaluations, expected_evaluations);
    }

    #[test]
    fn test_interpolates_poly() {
        let coeffs = vec![Fr::from(1), Fr::from(2), Fr::from(3), Fr::from(4)];

        let poly = UnivariatePoly::new(coeffs.clone());

        let evaluations = fft_evaluate(&poly);

        let interpolated_poly = fft_interpolate(&evaluations);

        assert_eq!(interpolated_poly.coefficient, coeffs);
    }
}
