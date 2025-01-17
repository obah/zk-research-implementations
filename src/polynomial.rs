use std::ops::{Add, Mul};

#[derive(Debug)]
struct UnivariatePoly {
    coefficient: Vec<f64>,
}

impl UnivariatePoly {
    fn new(coeff: Vec<f64>) -> Self {
        UnivariatePoly { coefficient: coeff }
    }

    fn evaluate(&self, x: u32) -> f64 {
        self.coefficient
            .iter()
            .enumerate()
            .map(|(index, coeff)| coeff * x.pow(index as u32) as f64)
            .sum()
    }

    fn degree(&self) -> usize {
        self.coefficient.len() - 1
    }

    fn scalar_mul(&self, scalar: f64) -> Self {
        let coefficients = self
            .coefficient
            .iter()
            .map(|point| point * scalar)
            .collect();

        UnivariatePoly::new(coefficients)
    }

    fn interpolate(points: Vec<(isize, isize)>) -> UnivariatePoly {
        let n = points.len();
        let mut result = UnivariatePoly::new(vec![0.0]);

        for i in 0..n {
            let (x_i, y_i) = points[i];
            let mut l_i = UnivariatePoly::new(vec![1.0]);

            for j in 0..n {
                if i != j {
                    let (x_j, _) = points[j];

                    let numerator = UnivariatePoly::new(vec![-x_j as f64, 1.0]);

                    let denominator = (x_i - x_j) as f64;

                    l_i = l_i * numerator.scalar_mul(1.0 / denominator);
                }
            }

            result = result + l_i.scalar_mul(y_i as f64);
        }

        result
    }
}

impl Add for UnivariatePoly {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let mut result = vec![0.0; self.coefficient.len().max(other.coefficient.len())];

        for (i, &coeff) in self.coefficient.iter().enumerate() {
            result[i] += coeff;
        }

        for (i, &coeff) in other.coefficient.iter().enumerate() {
            result[i] += coeff;
        }

        UnivariatePoly::new(result)
    }
}

impl Mul for UnivariatePoly {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        let mut coeffs = vec![0.0; self.degree() + other.degree() + 1];

        for (i, a) in self.coefficient.iter().enumerate() {
            for (j, b) in other.coefficient.iter().enumerate() {
                coeffs[i + j] += a * b;
            }
        }

        UnivariatePoly::new(coeffs)
    }
}

pub fn polynomial() {
    let points = vec![(0, 2), (1, 4), (2, 6)];
    let interpolated_poly = UnivariatePoly::interpolate(points);
    println!("Interpolated Polynomial: {:?}", interpolated_poly);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_returns_degree() {
        let poly_1 = UnivariatePoly {
            coefficient: vec![3.0, 4.0, 3.0],
        };

        assert!(poly_1.degree() == 2);
    }

    #[test]
    fn it_evaluates_poly() {
        let poly_1 = UnivariatePoly {
            coefficient: vec![3.0, 4.0, 3.0],
        };

        assert!(poly_1.evaluate(3) == 42.0);
    }

    #[test]
    fn it_scales() {
        let poly_1 = UnivariatePoly {
            coefficient: vec![3.0, 4.0, 3.0],
        };

        assert!(poly_1.scalar_mul(2.0).coefficient == vec![6.0, 8.0, 6.0]);
    }

    #[test]
    fn it_adds_correctly() {
        let poly_1 = UnivariatePoly {
            coefficient: vec![3.0, 4.0, 3.0],
        };

        let poly_2 = UnivariatePoly {
            coefficient: vec![-3.0, 0.0, 0.0, 4.0],
        };

        assert!((poly_1 + poly_2).coefficient == vec![0.0, 4.0, 3.0, 4.0]);
    }

    #[test]
    fn it_multiplies_correctly() {
        let poly_1 = UnivariatePoly {
            coefficient: vec![3.0, 4.0, 3.0],
        };

        let poly_2 = UnivariatePoly {
            coefficient: vec![-3.0, 0.0, 0.0, 4.0],
        };

        assert!((poly_1 * poly_2).coefficient == vec![-9.0, -12.0, -9.0, 12.0, 16.0, 12.0]);
    }

    #[test]
    fn it_interpolates_points() {
        let points = vec![(0, 2), (1, 4), (2, 6)];

        let new_poly = UnivariatePoly::interpolate(points);

        assert!(new_poly.coefficient == vec![2.0, 2.0, 0.0]);
    }
}
