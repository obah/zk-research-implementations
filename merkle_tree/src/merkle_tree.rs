use ark_ff::PrimeField;
use fiat_shamir::fiat_shamir_transcript::fq_vec_to_bytes;
use sha3::{Digest, Keccak256};

//new M.T
//insert leaf
//update leaf
//verify leaf existence
//proof leat existence
//hash

#[derive(PartialEq)]
enum LeafSide {
    Left,
    Right,
}

struct Leaf<F: PrimeField> {
    data_hash: F,
    leaf_id: Vec<usize>,
}

struct Node<F: PrimeField> {
    left_leaf: Leaf<F>,
    right_leaf: Leaf<F>,
    output_leaf: Leaf<F>,
}

struct MerkleTree<F: PrimeField> {
    nodes: Vec<Node<F>>,
    depth: usize,
}

impl<F: PrimeField> Leaf<F> {
    fn new(data: F, id: &[usize], is_hash: bool) -> Self {
        let leaf;

        if is_hash {
            leaf = Self {
                data_hash: data,
                leaf_id: id.to_vec(),
            }
        } else {
            let mut hasher = Keccak256::new();
            hasher.update(fq_vec_to_bytes(&[data]));
            let data_hash = F::from_le_bytes_mod_order(&hasher.finalize_reset());

            leaf = Self {
                data_hash,
                leaf_id: id.to_vec(),
            }
        }

        leaf
    }
}

impl<F: PrimeField> Node<F> {
    fn new(left_leaf: Leaf<F>, right_leaf: Leaf<F>) -> Self {
        assert_ne!(left_leaf.leaf_id, right_leaf.leaf_id, "invalid leaves");

        let mut hasher = Keccak256::new();

        hasher.update(fq_vec_to_bytes(&[left_leaf.data_hash]));
        hasher.update(fq_vec_to_bytes(&[right_leaf.data_hash]));
        let output_hash = F::from_le_bytes_mod_order(&hasher.finalize_reset());

        let mut new_id = left_leaf.leaf_id.clone();
        new_id.extend(&right_leaf.leaf_id);

        let output_leaf = Leaf::new(output_hash, &new_id, true);

        Self {
            left_leaf,
            right_leaf,
            output_leaf,
        }
    }

    fn update(&mut self, data_side: LeafSide, data: F) -> F {
        assert_eq!(self.left_leaf.leaf_id.len(), 1, "Can't update this node");

        let leaf_id_to_update = if data_side == LeafSide::Left {
            &self.left_leaf.leaf_id
        } else {
            &self.right_leaf.leaf_id
        };

        let new_leaf = Leaf::new(data, &leaf_id_to_update, false);

        if data_side == LeafSide::Left {
            self.left_leaf = new_leaf;
        } else {
            self.right_leaf = new_leaf;
        }

        let mut hasher = Keccak256::new();
        hasher.update(fq_vec_to_bytes(&[self.left_leaf.data_hash]));
        hasher.update(fq_vec_to_bytes(&[self.right_leaf.data_hash]));
        let new_output_hash = F::from_le_bytes_mod_order(&hasher.finalize_reset());

        let new_output_leaf = Leaf::new(new_output_hash, &self.output_leaf.leaf_id, true);

        self.output_leaf = new_output_leaf;

        self.output_leaf.data_hash
    }

    fn delete(&mut self, data_side: LeafSide) {
        self.update(data_side, F::zero());
    }
}

impl<F: PrimeField> MerkleTree<F> {
    fn new() -> Self {
        //todo add checks based on the depth and inputs
        todo!()
    }

    //?this should return a proof maybe as a path
    fn create_proof() {
        todo!()
    }

    //?this should take in the path
    fn verify() -> bool {
        todo!()
    }

    //? returns path also
    fn get_leaf_path() {
        todo!()
    }

    fn recompute_root_hash() {
        todo!()
    }
}

#[cfg(test)]
mod test {
    use ark_bn254::Fq;

    use super::*;

    #[test]
    fn test_creates_new_leaf() {
        let data = Fq::from(100);
        let leaf = Leaf::new(data, &[1], false);

        assert_eq!(leaf.leaf_id, &[1]);
    }

    #[test]
    fn test_creates_new_node() {
        let data = Fq::from(100);

        let leaf_1 = Leaf::new(data, &[1], false);
        let leaf_2 = Leaf::new(data, &[2], false);

        let node = Node::new(leaf_1, leaf_2);

        assert_eq!(node.output_leaf.leaf_id, &[1, 2]);
    }

    #[test]
    fn test_update_node() {
        let data = Fq::from(100);
        let leaf_1 = Leaf::new(data, &[1], false);
        let leaf_2 = Leaf::new(data, &[2], false);
        let mut node = Node::new(leaf_1, leaf_2);

        let initial_leaf_hash = node.left_leaf.data_hash;
        let initial_right_hash = node.right_leaf.data_hash;
        let initial_output_hash = node.output_leaf.data_hash;

        let new_data = Fq::from(200);
        node.update(LeafSide::Left, new_data);

        assert_ne!(initial_output_hash, node.output_leaf.data_hash);
        assert_ne!(initial_leaf_hash, node.left_leaf.data_hash);
        assert_eq!(initial_right_hash, node.right_leaf.data_hash);

        assert_eq!(node.output_leaf.leaf_id, &[1, 2]);
        assert_eq!(node.left_leaf.leaf_id, &[1]);
    }

    #[test]
    fn test_delete_leaf() {
        let data = Fq::from(100);
        let leaf_1 = Leaf::new(data, &[1], false);
        let leaf_2 = Leaf::new(data, &[2], false);
        let mut node = Node::new(leaf_1, leaf_2);
        let initial_right_hash = node.right_leaf.data_hash;
        let initial_output_hash = node.output_leaf.data_hash;

        node.delete(LeafSide::Right);

        assert_eq!(node.right_leaf.leaf_id, &[2]);
        assert_ne!(initial_output_hash, node.output_leaf.data_hash);
        assert_ne!(initial_right_hash, node.right_leaf.data_hash);
    }
}
