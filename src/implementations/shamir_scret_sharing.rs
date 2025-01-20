use crate::implementations::polynomial::UnivariatePoly;
use ark_bn254::Fq;
use ark_std::rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

fn recreate_polynomial(points: Vec<(Fq, Fq)>, threshold: usize) -> UnivariatePoly<Fq> {
    if points.len() < threshold {
        panic!("Not enough points to recreate polynomial");
    }

    let selected_points: Vec<(Fq, Fq)> = if points.len() > 3 {
        points[0..4].to_vec()
    } else {
        points.clone()
    };

    UnivariatePoly::interpolate(selected_points)
}

fn get_secret(poly: &UnivariatePoly<Fq>, x_point: Fq) -> Fq {
    let secret = poly.evaluate(Fq::from(x_point));

    secret
}

fn share_points(shares: usize, poly: &UnivariatePoly<Fq>) -> Vec<(Fq, Fq)> {
    let mut rng = StdRng::from_entropy();

    let mut shares: Vec<(Fq, Fq)> = vec![(Fq::from(0), Fq::from(0)); shares];

    for i in 0..shares.len() {
        let random_x_point: Fq = rng.gen();

        let y_point = poly.evaluate(random_x_point);

        shares[i] = (random_x_point, y_point);
    }

    shares
}

#[cfg(test)]
mod tests {
    use crate::implementations::polynomial::UnivariatePoly;
    use ark_bn254::Fq;

    use super::{get_secret, recreate_polynomial, share_points};

    #[test]
    fn it_recreates_polynomial_with_valid_points() {
        let points = vec![
            (Fq::from(1), Fq::from(-1)),
            (Fq::from(2), Fq::from(5)),
            (Fq::from(3), Fq::from(13)),
        ];
        let threshold = 3;

        let secret_poly = recreate_polynomial(points, threshold);

        assert_eq!(
            secret_poly.coefficient,
            vec![Fq::from(-5), Fq::from(3), Fq::from(1)]
        );
    }

    #[test]
    fn it_returns_right_secret() {
        let secret_poly = UnivariatePoly::new(vec![Fq::from(-5), Fq::from(3), Fq::from(1)]);

        let secret_data = get_secret(&secret_poly, Fq::from(0));

        assert_eq!(secret_data, Fq::from(-5));
    }

    #[test]
    fn it_share_points_correctly() {
        let secret_poly = UnivariatePoly::new(vec![Fq::from(-5), Fq::from(3), Fq::from(1)]);

        let shares = share_points(10, &secret_poly);

        assert_eq!(shares.len(), 10);
    }

    #[test]
    fn it_all_works_properly() {
        let secret_poly =
            UnivariatePoly::new(vec![Fq::from(-5), Fq::from(3), Fq::from(1), Fq::from(0)]);
        let secret_point = Fq::from(0);
        let secret_data = Fq::from(-5);

        let shares = share_points(10, &secret_poly);
        let collected_shares = shares[2..6].to_vec();
        let recreated_poly = recreate_polynomial(collected_shares, 3);
        let recovered_secret = get_secret(&recreated_poly, secret_point);

        assert_eq!(recreated_poly.coefficient, secret_poly.coefficient);
        assert_eq!(recovered_secret, secret_data);
    }
}
