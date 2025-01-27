use ark_bn254::Fq;
use ark_ff::PrimeField;
use criterion::{criterion_group, criterion_main, Criterion};
use rand::Rng;
use std::hint::black_box;
use std::ops::Range;

fn generate_random_field_elements<F: PrimeField>(range: Range<usize>) -> Vec<F> {
    let mut rng = rand::thread_rng();
    range.map(|_| F::from(rng.gen::<u64>())).collect()
}

fn bench_generate_bhc(c: &mut Criterion) {
    let evaluations = generate_random_field_elements::<Fr>(0..256);
    c.bench_function("generate_bhc", |b| {
        b.iter(|| BHC::generate_bhc(evaluations.clone()))
    });
}

fn bench_pair_points(c: &mut Criterion) {
    let evaluations = generate_random_field_elements::<Fr>(0..256);
    let bhc = BHC::generate_bhc(evaluations);
    c.bench_function("pair_points", |b| b.iter(|| bhc.pair_points(3)));
}

fn bench_partial_evaluate(c: &mut Criterion) {
    let evaluations = generate_random_field_elements::<Fr>(0..256);
    let poly = MultilinearPoly::new(evaluations);
    let value = Fr::from(5); // Example value
    c.bench_function("partial_evaluate", |b| {
        b.iter(|| poly.partial_evaluate(value, 3))
    });
}

fn bench_evaluate(c: &mut Criterion) {
    let evaluations = generate_random_field_elements::<Fr>(0..256);
    let poly = MultilinearPoly::new(evaluations);
    let values = generate_random_field_elements::<Fr>(0..8); // Example values
    c.bench_function("evaluate", |b| b.iter(|| poly.evaluate(values.clone())));
}

criterion_group!(
    benches,
    bench_generate_bhc,
    bench_pair_points,
    bench_partial_evaluate,
    bench_evaluate
);
criterion_main!(benches);
