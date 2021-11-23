
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rand::random;
use fast_kd::distance::squared_euclidean;
use fast_kd::KdTree;
use kdtree::KdTree as KdTree2;
use kdtree::distance::squared_euclidean as squared_euclidean2;

pub fn bench_insert_1M(c : &mut Criterion) {
    let mut tree : KdTree<f64, 2> = KdTree::new();
    for _ in 0..1_000_000 {
        tree.add(&[random(), random()]);
    }
    c.bench_function("tree add", |b|
        b.iter(||{
        tree.add(&[random(), random()])
    }));
}

// pub fn bench_insert_1M_kd(b : &mut Bencher) {
//     let mut tree : KdTree2<f64, i64, &[f64;2]> = KdTree2::new(2);
//     for _ in 0..1_000_000 {
//         tree.add(&[random(), random()], 0);
//     }
//     b.iter(|| {
//         tree.add(&[random(), random()], 0);
//     })
// }

pub fn bench_predict_1M_1000(c : &mut Criterion) {
    let mut tree : KdTree<f64, 2> = KdTree::new();
    for _ in 0..1_000_000 {
        tree.add(&[random(), random()]);
    }
    c.bench_function("tree predict", |b| b.iter(|| {
            tree.best_n_within(
                &[random(), random()],
                f64::INFINITY,
                1000,
                &squared_euclidean
            ).unwrap()
    }));
}

pub fn bench_predict_1M_1000_kd(c : &mut Criterion) {
    let mut tree : KdTree2<f64, i64, &[f64;2]> = KdTree2::new(2);
    let mut items:Vec<[f64;2]> = Vec::with_capacity(1_000_000);
    for _ in 0..1_000_000 {
        items.push([random(), random()]);
    }
    let items = items;
    for p in items.iter() {
        tree.add(p, 0);
    }
    c.bench_function("tree predict kd" , |b| b.iter(|| {
        tree.nearest(
            &[random(), random()],
            1000,
            &squared_euclidean2
        ).unwrap()
    }));
}


criterion_group!(
    benches,
    bench_insert_1M,
    //bench_insert_1M_kd,
    bench_predict_1M_1000,
    bench_predict_1M_1000_kd,
);
criterion_main!(benches);
