use ark_ff::PrimeField;
use fiat_shamir::fiat_shamir_transcript::fq_vec_to_bytes;
use sha3::{Digest, Keccak256};

#[derive(PartialEq, Copy, Clone, Debug)]
enum LeafSide {
    Left,
    Right,
}

struct LeafPath(Vec<Vec<usize>>);

#[derive(Clone, Debug)]
struct Leaf<F: PrimeField> {
    data_hash: F,
    leaf_id: Vec<usize>,
}

#[derive(Clone, Debug)]
struct Node<F: PrimeField> {
    left_leaf: Leaf<F>,
    right_leaf: Leaf<F>,
    output_leaf: Leaf<F>,
    node_id: usize,
    data_side: LeafSide,
}

#[derive(Debug)]
struct MerkleTree<F: PrimeField> {
    tree: Vec<Vec<Node<F>>>,
    depth: usize,
    next_leaf_id: usize,
}

#[derive(Clone)]
struct ProofData<F: PrimeField> {
    data_hash: F,
    data_side: LeafSide,
}

#[derive(Clone)]
struct MerkleProof<F: PrimeField> {
    data: F,
    proof: Vec<ProofData<F>>,
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
    fn new(left_leaf: Leaf<F>, right_leaf: Leaf<F>, node_id: usize, data_side: LeafSide) -> Self {
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
            node_id,
            data_side,
        }
    }

    fn update(&mut self, data_side: LeafSide, data: F, is_hash: bool) -> F {
        let leaf_id_to_update = if data_side == LeafSide::Left {
            &self.left_leaf.leaf_id
        } else {
            &self.right_leaf.leaf_id
        };

        let new_leaf = Leaf::new(data, &leaf_id_to_update, is_hash);

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
        self.update(data_side, F::zero(), true);
    }
}

impl<F: PrimeField> MerkleTree<F> {
    fn new(depth: usize) -> Self {
        let num_of_leaves = 1 << depth;

        let mut tree: Vec<Vec<Node<F>>> = Vec::with_capacity(depth);

        let mut leaves: Vec<Leaf<F>> = (0..num_of_leaves)
            .map(|index| Leaf::new(F::zero(), &[index], true))
            .collect();

        for _ in 0..depth {
            let mut nodes = Vec::new();

            for (idx, pair) in leaves.chunks(2).enumerate() {
                let data_side = if idx % 2 == 1 {
                    LeafSide::Right
                } else {
                    LeafSide::Left
                };

                let node = Node::new(pair[0].clone(), pair[1].clone(), idx, data_side);
                nodes.push(node);
            }

            tree.push(nodes.clone());
            leaves = nodes.iter().map(|node| node.output_leaf.clone()).collect();
        }

        Self {
            tree,
            depth,
            next_leaf_id: 0,
        }
    }

    fn insert_leaf(&mut self, data: F, is_hash: bool) {
        let id = self.next_leaf_id;

        if id == 1 << self.depth {
            panic!("All leaves are filled!");
        }

        let new_leaf = Leaf::new(data, &[id], is_hash);

        self.tree[0].iter_mut().for_each(|node| {
            if node.left_leaf.leaf_id[0] == id {
                node.left_leaf = new_leaf.clone();
                return;
            } else if node.right_leaf.leaf_id[0] == id {
                node.right_leaf = new_leaf.clone();
                return;
            }
        });

        self.next_leaf_id += 1;

        self.recompute_root_hash(id);
    }

    fn update_leaf(&mut self, data: F, leaf_id: usize, is_hash: bool) {
        assert!(leaf_id < 1 << self.depth, "Invalid leaf id");

        let new_leaf = Leaf::new(data, &[leaf_id], is_hash);

        self.tree[0].iter_mut().for_each(|node| {
            if node.left_leaf.leaf_id[0] == leaf_id {
                node.left_leaf = new_leaf.clone();
                return;
            } else if node.right_leaf.leaf_id[0] == leaf_id {
                node.right_leaf = new_leaf.clone();
                return;
            }
        });

        self.recompute_root_hash(leaf_id);
    }

    fn get_root_hash(&self) -> F {
        let last_node = self.depth - 1;

        self.tree[last_node][0].output_leaf.data_hash
    }

    fn create_proof(&self, data_to_prove: F, data_id: usize) -> MerkleProof<F> {
        let selected_node = self.tree[0].iter().find(|node| {
            node.left_leaf.leaf_id[0] == data_id || node.right_leaf.leaf_id[0] == data_id
        });

        if selected_node.is_none() {
            panic!("Data with given id not found");
        }

        let mut selected_node = selected_node.unwrap();

        let selected_leaf;
        let proof_leaf;
        let proof_side;

        if selected_node.left_leaf.leaf_id[0] == data_id {
            selected_leaf = selected_node.left_leaf.clone();
            proof_leaf = selected_node.right_leaf.clone();
            proof_side = LeafSide::Right;
        } else {
            selected_leaf = selected_node.right_leaf.clone();
            proof_leaf = selected_node.left_leaf.clone();
            proof_side = LeafSide::Left;
        }

        let mut hasher = Keccak256::new();
        hasher.update(fq_vec_to_bytes(&[data_to_prove]));
        let data_to_prove_hash = F::from_le_bytes_mod_order(&hasher.finalize_reset());

        assert_eq!(
            selected_leaf.data_hash, data_to_prove_hash,
            "invalid leaf data"
        );

        let mut proofs: Vec<ProofData<F>> = Vec::with_capacity(self.depth + 1);

        let first_proof = ProofData {
            data_hash: proof_leaf.data_hash,
            data_side: proof_side,
        };

        proofs.push(first_proof);

        let mut current_side = selected_node.data_side;
        for idx in 0..(self.depth - 1) {
            let parent_node = self.tree[idx + 1]
                .iter()
                .find(|node| {
                    node.left_leaf.leaf_id == selected_node.output_leaf.leaf_id
                        || node.right_leaf.leaf_id == selected_node.output_leaf.leaf_id
                })
                .unwrap();

            let data_side;
            let proof_hash;

            if current_side == LeafSide::Left {
                proof_hash = parent_node.right_leaf.data_hash;
                data_side = LeafSide::Right;
            } else {
                proof_hash = parent_node.left_leaf.data_hash;
                data_side = LeafSide::Left;
            }

            let proof = ProofData {
                data_hash: proof_hash,
                data_side,
            };

            proofs.push(proof);

            current_side = parent_node.data_side;
            selected_node = parent_node;
        }

        MerkleProof {
            data: data_to_prove,
            proof: proofs,
        }
    }

    fn verify(&self, proof: MerkleProof<F>) -> bool {
        let root_hash = self.get_root_hash();

        let mut hasher = Keccak256::new();

        hasher.update(fq_vec_to_bytes(&[proof.data]));
        let data_to_prove_hash = F::from_le_bytes_mod_order(&hasher.finalize_reset());

        let mut current_hash = data_to_prove_hash;
        let mut verified = false;

        for proof_data in proof.proof {
            if proof_data.data_side == LeafSide::Left {
                hasher.update(fq_vec_to_bytes(&[proof_data.data_hash]));
                hasher.update(fq_vec_to_bytes(&[current_hash]));
            } else {
                hasher.update(fq_vec_to_bytes(&[current_hash]));
                hasher.update(fq_vec_to_bytes(&[proof_data.data_hash]));
            }

            current_hash = F::from_le_bytes_mod_order(&hasher.finalize_reset());
        }

        if root_hash == current_hash {
            verified = true;
        }

        verified
    }

    fn get_leaf_path(&self, leaf_id: usize) -> LeafPath {
        let all_ids: Vec<Vec<Vec<usize>>> = self
            .tree
            .iter()
            .map(|layer| {
                layer
                    .iter()
                    .map(|node| node.output_leaf.leaf_id.clone())
                    .collect()
            })
            .collect();

        let affected_path = all_ids
            .iter()
            .flatten()
            .filter(|id| id.contains(&leaf_id))
            .cloned()
            .collect();

        LeafPath(affected_path)
    }

    fn recompute_root_hash(&mut self, mut affected_leaf_id: usize) {
        assert!(affected_leaf_id < 1 << self.depth, "Invalid leaf id");

        let affected_leaf_path = self.get_leaf_path(affected_leaf_id);

        for (index, id) in affected_leaf_path.0.iter().enumerate() {
            let mut affected_node: Node<F> = self.tree[index]
                .iter()
                .find(|node| node.output_leaf.leaf_id == *id)
                .unwrap()
                .clone();

            let data_side;
            let affected_data;

            if affected_leaf_id % 2 == 1 {
                data_side = LeafSide::Right;
                affected_data = affected_node.right_leaf.data_hash;
            } else {
                data_side = LeafSide::Left;
                affected_data = affected_node.left_leaf.data_hash;
            }

            let new_hash = affected_node.update(data_side, affected_data, true);

            affected_leaf_id = affected_node.node_id;

            self.tree[index].iter_mut().for_each(|node| {
                if node.output_leaf.leaf_id == *id {
                    node.output_leaf.data_hash = new_hash;
                }
            });

            let next_leaf_id = affected_node.output_leaf.leaf_id;

            if index < self.depth - 1 {
                self.tree[index + 1].iter_mut().for_each(|node| {
                    if node.left_leaf.leaf_id == next_leaf_id {
                        node.left_leaf.data_hash = new_hash;
                    } else {
                        node.right_leaf.data_hash = new_hash;
                    }
                })
            }
        }
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

        let node = Node::new(leaf_1, leaf_2, 0, LeafSide::Left);

        assert_eq!(node.output_leaf.leaf_id, &[1, 2]);
    }

    #[test]
    fn test_update_node() {
        let data = Fq::from(100);
        let leaf_1 = Leaf::new(data, &[1], false);
        let leaf_2 = Leaf::new(data, &[2], false);
        let mut node = Node::new(leaf_1, leaf_2, 0, LeafSide::Left);

        let initial_leaf_hash = node.left_leaf.data_hash;
        let initial_right_hash = node.right_leaf.data_hash;
        let initial_output_hash = node.output_leaf.data_hash;

        let new_data = Fq::from(200);
        node.update(LeafSide::Left, new_data, false);

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
        let mut node = Node::new(leaf_1, leaf_2, 0, LeafSide::Left);
        let initial_right_hash = node.right_leaf.data_hash;
        let initial_output_hash = node.output_leaf.data_hash;

        node.delete(LeafSide::Right);

        assert_eq!(node.right_leaf.data_hash, Fq::from(0));
        assert_eq!(node.right_leaf.leaf_id, &[2]);
        assert_ne!(initial_output_hash, node.output_leaf.data_hash);
        assert_ne!(initial_right_hash, node.right_leaf.data_hash);
    }

    #[test]
    fn test_create_tree() {
        let depth = 2;
        let merkle_tree: MerkleTree<Fq> = MerkleTree::new(depth);

        assert_eq!(merkle_tree.tree.len(), depth);
        assert_eq!(merkle_tree.tree[0][0].left_leaf.data_hash, Fq::from(0));
    }

    #[test]
    fn test_get_root_hash() {
        let depth = 2;
        let merkle_tree: MerkleTree<Fq> = MerkleTree::new(depth);
        let root_leaf = merkle_tree.tree[depth - 1][0].output_leaf.clone();

        assert_eq!(root_leaf.leaf_id, &[0, 1, 2, 3]);
        assert_eq!(merkle_tree.tree[depth - 1].len(), 1);

        let mut hasher = Keccak256::new();
        hasher.update(fq_vec_to_bytes(&[Fq::from(0)]));
        hasher.update(fq_vec_to_bytes(&[Fq::from(0)]));
        let hash_1 = Fq::from_le_bytes_mod_order(&hasher.finalize_reset());

        hasher.update(fq_vec_to_bytes(&[hash_1]));
        hasher.update(fq_vec_to_bytes(&[hash_1]));
        let root_hash = Fq::from_le_bytes_mod_order(&hasher.finalize_reset());

        assert_eq!(root_hash, merkle_tree.get_root_hash());
    }

    #[test]
    fn test_get_leaf_path() {
        let depth = 2;
        let merkle_tree: MerkleTree<Fq> = MerkleTree::new(depth);
        let leaf_path = merkle_tree.get_leaf_path(0);

        let expected_affected_path = vec![vec![0, 1], vec![0, 1, 2, 3]];

        assert_eq!(leaf_path.0, expected_affected_path);
    }

    #[test]
    fn test_recompute_root_hash() {
        let depth = 2;
        let mut merkle_tree: MerkleTree<Fq> = MerkleTree::new(depth);

        let new_data = Fq::from(10);
        let new_leaf = Leaf::new(new_data, &[0], false);

        merkle_tree.tree[0][0].left_leaf = new_leaf;

        merkle_tree.recompute_root_hash(0);

        let mut hasher = Keccak256::new();
        hasher.update(fq_vec_to_bytes(&[new_data]));
        let new_data_hash = Fq::from_le_bytes_mod_order(&hasher.finalize_reset());

        hasher.update(fq_vec_to_bytes(&[new_data_hash]));
        hasher.update(fq_vec_to_bytes(&[Fq::from(0)]));
        let hash_1 = Fq::from_le_bytes_mod_order(&hasher.finalize_reset());

        hasher.update(fq_vec_to_bytes(&[Fq::from(0)]));
        hasher.update(fq_vec_to_bytes(&[Fq::from(0)]));
        let hash_2 = Fq::from_le_bytes_mod_order(&hasher.finalize_reset());

        hasher.update(fq_vec_to_bytes(&[hash_1]));
        hasher.update(fq_vec_to_bytes(&[hash_2]));
        let expected_root_hash = Fq::from_le_bytes_mod_order(&hasher.finalize_reset());

        assert_eq!(expected_root_hash, merkle_tree.get_root_hash());
    }

    #[test]
    fn test_insert_leaf() {
        let depth = 2;
        let mut merkle_tree: MerkleTree<Fq> = MerkleTree::new(depth);

        let new_data = Fq::from(10);
        merkle_tree.insert_leaf(new_data, false);

        let mut hasher = Keccak256::new();
        hasher.update(fq_vec_to_bytes(&[new_data]));
        let new_data_hash = Fq::from_le_bytes_mod_order(&hasher.finalize_reset());

        assert_eq!(new_data_hash, merkle_tree.tree[0][0].left_leaf.data_hash);

        hasher.update(fq_vec_to_bytes(&[new_data_hash]));
        hasher.update(fq_vec_to_bytes(&[Fq::from(0)]));
        let hash_1 = Fq::from_le_bytes_mod_order(&hasher.finalize_reset());

        hasher.update(fq_vec_to_bytes(&[Fq::from(0)]));
        hasher.update(fq_vec_to_bytes(&[Fq::from(0)]));
        let hash_2 = Fq::from_le_bytes_mod_order(&hasher.finalize_reset());

        hasher.update(fq_vec_to_bytes(&[hash_1]));
        hasher.update(fq_vec_to_bytes(&[hash_2]));
        let expected_root_hash = Fq::from_le_bytes_mod_order(&hasher.finalize_reset());

        assert_eq!(expected_root_hash, merkle_tree.get_root_hash());
    }

    #[test]
    fn test_update_leaf() {
        let depth = 2;
        let mut merkle_tree: MerkleTree<Fq> = MerkleTree::new(depth);

        let new_data = Fq::from(10);
        merkle_tree.update_leaf(new_data, 1, false);

        let mut hasher = Keccak256::new();
        hasher.update(fq_vec_to_bytes(&[new_data]));
        let new_data_hash = Fq::from_le_bytes_mod_order(&hasher.finalize_reset());

        assert_eq!(new_data_hash, merkle_tree.tree[0][0].right_leaf.data_hash);

        hasher.update(fq_vec_to_bytes(&[Fq::from(0)]));
        hasher.update(fq_vec_to_bytes(&[new_data_hash]));
        let hash_1 = Fq::from_le_bytes_mod_order(&hasher.finalize_reset());

        hasher.update(fq_vec_to_bytes(&[Fq::from(0)]));
        hasher.update(fq_vec_to_bytes(&[Fq::from(0)]));
        let hash_2 = Fq::from_le_bytes_mod_order(&hasher.finalize_reset());

        hasher.update(fq_vec_to_bytes(&[hash_1]));
        hasher.update(fq_vec_to_bytes(&[hash_2]));
        let expected_root_hash = Fq::from_le_bytes_mod_order(&hasher.finalize_reset());

        assert_eq!(expected_root_hash, merkle_tree.get_root_hash());
    }

    #[test]
    fn test_proof_and_verify() {
        let depth = 3;
        let mut merkle_tree: MerkleTree<Fq> = MerkleTree::new(depth);

        let new_data = Fq::from(10);
        merkle_tree.insert_leaf(new_data, false);

        let proof = merkle_tree.create_proof(new_data, 0);

        let is_verified = merkle_tree.verify(proof);

        assert_eq!(is_verified, true);
    }

    #[test]
    fn test_verify_invalid_proof() {
        let depth = 2;
        let merkle_tree: MerkleTree<Fq> = MerkleTree::new(depth);
        let new_data = Fq::from(10);

        let fake_proof_data = ProofData {
            data_hash: Fq::from(0),
            data_side: LeafSide::Left,
        };

        let invalid_proof = MerkleProof {
            data: new_data,
            proof: vec![fake_proof_data; 2],
        };

        let is_verified = merkle_tree.verify(invalid_proof);
        assert_eq!(is_verified, false);
    }
}
