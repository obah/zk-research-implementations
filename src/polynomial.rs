use std::ops::{Add, Mul};

#[derive(Debug)]
struct UnivariatePoly {
    coefficient: Vec<f64>,
}

struct Point(isize, isize);

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

    fn scalar_mul(&self, val: f64) -> Self {
        let coefficients = self.coefficient.iter().map(|point| point * val).collect();

        UnivariatePoly::new(coefficients)
    }

    fn interpolate(points: Vec<Point>) -> UnivariatePoly {
        todo!()
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
}
