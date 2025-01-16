struct UnivariatePoly {
    coefficient: Vec<usize>,
}

// struct Point {
//     x: usize,
//     y: usize,
// }

impl UnivariatePoly {
    fn evaluate(&self, x: usize) -> usize {
        self.coefficient
            .iter()
            .enumerate()
            .map(|(index, coeff)| coeff * x.pow(index.try_into().unwrap()))
            .sum()

        //map each item to coeff * var ^^ index then sum
    }

    fn degree(&self) -> usize {
        self.coefficient.len() - 1
    }

    // fn interpolate(points: Vec<Point>) -> UnivariatePoly {}
}

struct SparseUnivariatePoly {
    coefficient: Vec<(usize, usize)>,
}

impl SparseUnivariatePoly {
    fn evaluate(&self, x: usize) -> usize {
        self.coefficient
            .iter()
            .map(|coeff| coeff.0 * x.pow(coeff.1.try_into().unwrap()))
            .sum()
    }

    fn degree(&self) -> usize {
        let coefficients = &self.coefficient;

        let mut largest = coefficients[0].1;

        for (_, coeff) in coefficients.iter().enumerate() {
            let current_degree = coeff.1;

            if largest < current_degree {
                largest = current_degree;
            }
        }

        largest
    }
}

// {}

pub fn polynomial() {
    let polynomial1 = UnivariatePoly {
        coefficient: vec![1, 2, 4],
    };

    let result1 = polynomial1.evaluate(3);
    let degree1 = polynomial1.degree();

    println!("evaluation of polynomial1 is {result1} and its degree is {degree1}");

    let polynomial2 = SparseUnivariatePoly {
        coefficient: vec![(1, 0), (2, 1), (4, 2)],
    };

    let result2 = polynomial2.evaluate(3);
    let degree2 = polynomial2.degree();

    println!("evaluation of polynomial2 is {result2} and its degree is {degree2}");
}
