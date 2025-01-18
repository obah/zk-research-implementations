use crate::implementations::polynomial::UnivariatePoly;
use rand::Rng;

fn recreate_polynomial(points: Vec<(f64, f64)>, threshold: usize) -> UnivariatePoly {
    if points.len() < threshold {
        panic!("Not enough points to recreate polynomial");
    }

    let selected_points = if points.len() > 3 {
        points[0..4].to_vec()
    } else {
        points.clone()
    };

    UnivariatePoly::interpolate(selected_points)
}

fn get_secret(poly: &UnivariatePoly, x_point: u32) -> f64 {
    let secret = poly.evaluate(x_point);

    secret
}

fn share_points(shares: usize, poly: &UnivariatePoly) -> Vec<(f64, f64)> {
    let mut rng = rand::thread_rng();

    let mut shares: Vec<(f64, f64)> = vec![(0.0, 0.0); shares];

    for i in 0..shares.len() {
        let random_x_point: i32 = rng.gen_range(0..=100);

        let y_point = poly.evaluate(random_x_point as u32);

        shares[i] = (random_x_point as f64, y_point);
    }

    shares
}

#[cfg(test)]
mod tests {
    use crate::implementations::polynomial::UnivariatePoly;

    use super::{get_secret, recreate_polynomial, share_points};

    #[test]
    fn it_recreates_polynomial_with_valid_points() {
        let points = vec![(1.0, -1.0), (2.0, 5.0), (3.0, 13.0)];
        let threshold = 3;

        let secret_poly = recreate_polynomial(points, threshold);

        assert_eq!(secret_poly.coefficient, vec![-5.0, 3.0, 1.0]);
    }

    #[test]
    fn it_returns_right_secret() {
        let secret_poly = UnivariatePoly::new(vec![-5.0, 3.0, 1.0]);

        let secret_data = get_secret(&secret_poly, 0);

        assert_eq!(secret_data, -5.0);
    }

    #[test]
    fn it_share_points_correctly() {
        let secret_poly = UnivariatePoly::new(vec![-5.0, 3.0, 1.0]);

        let shares = share_points(10, &secret_poly);

        assert_eq!(shares.len(), 10);
    }

    #[test]
    fn it_all_works_properly() {
        let secret_poly = UnivariatePoly::new(vec![-5.0, 3.0, 1.0]);
        let secret_point = 0;
        let secret_data = -5.0;

        let shares = share_points(10, &secret_poly);
        let recreated_poly = recreate_polynomial(shares, 3);
        let recovered_secret = get_secret(&recreated_poly, secret_point);

        assert_eq!(recreated_poly.coefficient, secret_poly.coefficient);
        assert_eq!(recovered_secret.round(), secret_data);
    }
}
