use std::collections::{BinaryHeap};
use std::hash::Hash;

use num_traits::{Float, One, Zero};


use crate::heap_element::HeapElement;
use crate::util::distance_to_space;


#[derive(Clone, Debug)]
pub struct KdTree<A, const K: usize> {
    size: usize,


    min_bounds: [A; K],
    max_bounds: [A; K],
    content: Node<A, K>,
}

#[derive(Clone, Debug)]
pub enum Node<A, const K: usize> {
    Stem {
        left: Box<KdTree<A, K>>,
        right: Box<KdTree<A, K>>,
        split_value: A,
        split_dimension: usize,
    },
    Leaf {
        //points: Vec<[A; K]>,
        //bucket: Vec<T>,
        bucket: Vec<[A; K]>,
        capacity: usize,
    },
}

#[derive(Debug, PartialEq)]
pub enum ErrorKind {
    NonFiniteCoordinate,
    ZeroCapacity,
    Empty,
}

impl<A: Float + Zero + One, const K: usize> KdTree<A, K> {
    pub fn new() -> Self {
        KdTree::with_per_node_capacity(16).unwrap()
    }

    pub fn with_per_node_capacity(capacity: usize) -> Result<Self, ErrorKind> {
        if capacity == 0 {
            return Err(ErrorKind::ZeroCapacity);
        }

        Ok(KdTree {
            size: 0,
            min_bounds: [A::infinity(); K],
            max_bounds: [A::neg_infinity(); K],
            content: Node::Leaf {
                bucket: Vec::with_capacity(capacity + 1),
                capacity,
            },
        })
    }

    pub fn size(&self) -> usize {
        self.size
    }

    pub fn is_leaf(&self) -> bool {
        match &self.content {
            Node::Leaf { .. } => true,
            Node::Stem { .. } => false,
        }
    }


    pub fn best_n_within<F>(
        &self,
        point: &[A; K],
        radius: A,
        max_qty: usize,
        distance: &F,
    ) -> Result<Vec<&[A;K]>, ErrorKind>
    where
        F: Fn(&[A; K], &[A; K]) -> A,
    {
        if self.size == 0 {
            return Ok(vec![]);
        }

        self.check_point(point)?;

        let mut pending = BinaryHeap::new();
        let mut evaluated = BinaryHeap::<HeapElement<A, &[A;K]>>::with_capacity(self.size().min(max_qty + 1));
        let mut max_ev_dist = A::infinity();
        pending.push(HeapElement {
            distance: A::zero(),
            element: self,
        });

        while !pending.is_empty() {
            let curr = pending.pop().unwrap();
            if evaluated.len() == max_qty && -curr.distance > max_ev_dist {
                break;
            }
            let curr = curr.element;
            match curr.content {
                Node::Leaf {
                    ref bucket,
                    ..
                } => {
                    for p in bucket.iter() {
                        let d : A = distance(point, p);
                        let heap_elem = HeapElement {
                            distance: d,
                            element: p,
                        };

                        if evaluated.len() < max_qty {
                            evaluated.push(heap_elem);
                            max_ev_dist = evaluated.peek().unwrap().distance;
                        } else if max_ev_dist > heap_elem.distance {
                            evaluated.push(heap_elem);
                            evaluated.pop();
                            max_ev_dist = evaluated.peek().unwrap().distance;
                        }
                    }
                }
                Node::Stem {
                    ref left,
                    ref right,
                    ..
                } => {
                    let d_left :A = distance_to_space(
                        point,
                        &left.min_bounds,
                        &left.max_bounds,
                        distance
                    );
                    if d_left < radius {
                        pending.push(HeapElement {
                            distance: -d_left,
                            element: left,
                        });
                    }
                    let d_right:A = distance_to_space(
                        point,
                        &right.min_bounds,
                        &right.max_bounds,
                        distance
                    );
                    if d_right < radius {
                        pending.push(HeapElement {
                            distance: -d_right,
                            element: right,
                        });
                    }
                }
            }
        }

        Ok(evaluated.iter().map(|e| e.element).collect())
    }


    /// Add an element to the tree. The first argument specifies the location in kd space
    /// at which the element is located. The second argument is the data associated with
    /// that point in space.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use kiddo::KdTree;
    ///
    /// let mut tree: KdTree<f64, 3> = KdTree::new();
    ///
    /// tree.add(&[1.0, 2.0, 5.0])?;
    /// tree.add(&[1.1, 2.1, 5.1])?;
    ///
    /// assert_eq!(tree.size(), 2);
    /// ```
    pub fn add(&mut self, point: &[A; K]) -> Result<(), ErrorKind> {
        self.check_point(point)?;
        self.add_unchecked(point);
        Ok(())
    }

    fn add_unchecked(&mut self, point: &[A; K]) {
        let res = match &mut self.content {
            Node::Leaf { .. } => {
                self.add_to_bucket(point);
            }

            Node::Stem {
                left,
                right,
                split_dimension,
                split_value,
            } => {
                if point[*split_dimension] < *split_value {
                    // belongs_in_left
                    left.add_unchecked(point)
                } else {
                    right.add_unchecked(point)
                }
            }
        };

        self.extend(point);
        self.size += 1;
    }

    fn add_to_bucket(&mut self, point: &[A; K]) {
        self.extend(point);
        let cap;
        match &mut self.content {
            Node::Leaf {
                bucket,
                capacity,
            } => {
                bucket.push(*point);
                cap = *capacity;
            }
            Node::Stem { .. } => unreachable!(),
        }

        self.size += 1;
        if self.size > cap {
            self.split();
        }
    }

    fn split(&mut self) {
        match &mut self.content {
            Node::Leaf {
                bucket,
                capacity,
            } => {
                let mut split_dimension:usize = 0;
                let mut max = A::zero();
                for dim in 0..K {
                    let diff = self.max_bounds[dim] - self.min_bounds[dim];
                    if !diff.is_nan() && diff > max {
                        max = diff;
                        split_dimension = dim;
                    }
                }

                let split_value = self.min_bounds[split_dimension] + max / A::from(2.0).unwrap();
                let mut left = Box::new(KdTree::with_per_node_capacity(*capacity).unwrap());
                let mut right = Box::new(KdTree::with_per_node_capacity(*capacity).unwrap());

                while !bucket.is_empty() {
                    let point= bucket.pop().unwrap();
                    if point[split_dimension] < split_value {
                        // belongs_in_left
                        left.add_to_bucket(&point);
                    } else {
                        right.add_to_bucket(&point);
                    }
                }

                self.content = Node::Stem {
                    left,
                    right,
                    split_value,
                    split_dimension,
                }
            }
            Node::Stem { .. } => unreachable!(),
        }
    }

    fn extend(&mut self, point: &[A; K]) {
        let min = self.min_bounds.iter_mut();
        let max = self.max_bounds.iter_mut();
        for ((l, h), v) in min.zip(max).zip(point.iter()) {
            if v < l {
                *l = *v
            }
            if v > h {
                *h = *v
            }
        }
    }

    fn check_point(&self, point: &[A; K]) -> Result<(), ErrorKind> {
        if point.iter().all(|n| n.is_finite()) {
            Ok(())
        } else {
            Err(ErrorKind::NonFiniteCoordinate)
        }
    }
}


impl std::fmt::Display for ErrorKind {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "KdTree error: {}", self)
    }
}