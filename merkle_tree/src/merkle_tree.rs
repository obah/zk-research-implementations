use ark_ff::PrimeField;
use fiat_shamir::fiat_shamir_transcript::fq_vec_to_bytes;
use sha3::{Digest, Keccak256};
use std::error::Error;

#[derive(PartialEq, Copy, Clone, Debug)]
enum LeafSide {
    Left,
    Right,
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

#[derive(Debug)]
struct MerkleTree<F: PrimeField> {
    leaves: Vec<F>,
    tree: Vec<Vec<F>>,
    depth: usize,
}

impl<F: PrimeField> MerkleTree<F> {
    pub fn new(depth: usize) -> Self {
        let num_leaves = 1 << depth;
        let leaves = vec![F::zero(); num_leaves];
        let mut tree = Vec::with_capacity(depth);
        let mut current_level = leaves.clone();

        for _ in 0..depth {
            let next_level = current_level
                .chunks(2)
                .map(|pair| Self::hash_pair(pair[0], pair[1]))
                .collect::<Vec<_>>();
            tree.push(next_level.clone());
            current_level = next_level;
        }

        Self {
            leaves,
            tree,
            depth,
        }
    }

    pub fn new_with_inputs(depth: usize, inputs: Vec<F>) -> Result<Self, Box<dyn Error>> {
        let num_leaves = 1 << depth;
        if inputs.len() > num_leaves {
            return Err("Too many inputs for tree depth".into());
        }

        let mut leaves = vec![F::zero(); num_leaves];
        for (i, input) in inputs.iter().enumerate() {
            leaves[i] = Self::compute_hash(*input);
        }

        let mut tree = Vec::with_capacity(depth);
        let mut current_level = leaves.clone();

        for _ in 0..depth {
            let next_level = current_level
                .chunks(2)
                .map(|pair| Self::hash_pair(pair[0], pair[1]))
                .collect::<Vec<_>>();
            tree.push(next_level.clone());
            current_level = next_level;
        }

        Ok(Self {
            leaves,
            tree,
            depth,
        })
    }

    pub fn update_leaf(
        &mut self,
        leaf_id: usize,
        data: F,
        is_hash: bool,
    ) -> Result<(), Box<dyn Error>> {
        if leaf_id >= 1 << self.depth {
            return Err("Invalid leaf ID".into());
        }

        let new_hash = if is_hash {
            data
        } else {
            Self::compute_hash(data)
        };

        self.leaves[leaf_id] = new_hash;
        self.recompute_path(leaf_id);

        Ok(())
    }

    fn recompute_path(&mut self, leaf_id: usize) {
        let mut current_hash = self.leaves[leaf_id];
        let mut index = leaf_id;

        for level in 0..self.depth {
            let sibling_index = index ^ 1;

            let sibling_hash = if level == 0 {
                self.leaves[sibling_index]
            } else {
                self.tree[level - 1][sibling_index]
            };

            let (left, right) = if index % 2 == 0 {
                (current_hash, sibling_hash)
            } else {
                (sibling_hash, current_hash)
            };

            current_hash = Self::hash_pair(left, right);

            let parent_index = index / 2;

            self.tree[level][parent_index] = current_hash;
            index = parent_index;
        }
    }

    pub fn get_root_hash(&self) -> F {
        self.tree[self.depth - 1][0]
    }

    pub fn create_proof(
        &self,
        data_to_prove: F,
        leaf_id: usize,
    ) -> Result<MerkleProof<F>, Box<dyn Error>> {
        if leaf_id >= 1 << self.depth {
            return Err("Invalid leaf ID".into());
        }

        let data_hash = Self::compute_hash(data_to_prove);

        if self.leaves[leaf_id] != data_hash {
            return Err("Data does not match the leaf hash".into());
        }

        let mut proof = Vec::with_capacity(self.depth);
        let mut index = leaf_id;

        for level in 0..self.depth {
            let sibling_index = index ^ 1;

            let sibling_hash = if level == 0 {
                self.leaves[sibling_index]
            } else {
                self.tree[level - 1][sibling_index]
            };

            let data_side = if index % 2 == 0 {
                LeafSide::Right
            } else {
                LeafSide::Left
            };

            proof.push(ProofData {
                data_hash: sibling_hash,
                data_side,
            });

            index /= 2;
        }

        Ok(MerkleProof {
            data: data_to_prove,
            proof,
        })
    }

    pub fn verify(&self, proof: MerkleProof<F>) -> bool {
        let root_hash = self.get_root_hash();
        let mut current_hash = Self::compute_hash(proof.data);

        for proof_data in proof.proof {
            let (left, right) = match proof_data.data_side {
                LeafSide::Left => (proof_data.data_hash, current_hash),
                LeafSide::Right => (current_hash, proof_data.data_hash),
            };

            current_hash = Self::hash_pair(left, right);
        }

        root_hash == current_hash
    }

    fn compute_hash(data: F) -> F {
        let mut hasher = Keccak256::new();
        hasher.update(fq_vec_to_bytes(&[data]));

        F::from_le_bytes_mod_order(&hasher.finalize_reset())
    }

    fn hash_pair(left: F, right: F) -> F {
        let mut hasher = Keccak256::new();
        hasher.update(fq_vec_to_bytes(&[left]));
        hasher.update(fq_vec_to_bytes(&[right]));

        F::from_le_bytes_mod_order(&hasher.finalize_reset())
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bn254::Fq;

    #[test]
    fn test_create_tree() {
        let depth = 2;
        let merkle_tree: MerkleTree<Fq> = MerkleTree::new(depth);
        assert_eq!(merkle_tree.leaves.len(), 4);
        assert_eq!(merkle_tree.tree.len(), depth);
        assert_eq!(merkle_tree.tree[0].len(), 2);
        assert_eq!(merkle_tree.tree[1].len(), 1);

        for leaf in merkle_tree.leaves.iter() {
            assert_eq!(*leaf, Fq::from(0));
        }

        let zero_hash = Fq::from(0);
        let hash_1 = MerkleTree::<Fq>::hash_pair(zero_hash, zero_hash);
        let hash_2 = MerkleTree::<Fq>::hash_pair(hash_1, hash_1);

        assert_eq!(merkle_tree.get_root_hash(), hash_2);
    }

    #[test]
    fn test_get_root_hash() {
        let depth = 2;
        let merkle_tree: MerkleTree<Fq> = MerkleTree::new(depth);
        let zero_hash = Fq::from(0);
        let hash_1 = MerkleTree::<Fq>::hash_pair(zero_hash, zero_hash);
        let hash_2 = MerkleTree::<Fq>::hash_pair(hash_1, hash_1);

        assert_eq!(merkle_tree.get_root_hash(), hash_2);
    }

    #[test]
    fn test_update_leaf() {
        let depth = 2;
        let mut merkle_tree: MerkleTree<Fq> = MerkleTree::new(depth);
        let new_data = Fq::from(10);

        merkle_tree.update_leaf(1, new_data, false).unwrap();

        let hash_new_data = MerkleTree::<Fq>::compute_hash(new_data);

        assert_eq!(merkle_tree.leaves[1], hash_new_data);

        let hash_0 = Fq::from(0);
        let hash_2 = Fq::from(0);
        let hash_3 = Fq::from(0);
        let hash_01 = MerkleTree::<Fq>::hash_pair(hash_0, hash_new_data);
        let hash_23 = MerkleTree::<Fq>::hash_pair(hash_2, hash_3);
        let expected_root_hash = MerkleTree::<Fq>::hash_pair(hash_01, hash_23);

        assert_eq!(merkle_tree.get_root_hash(), expected_root_hash);
    }

    #[test]
    fn test_delete_leaf() {
        let depth = 2;
        let mut merkle_tree: MerkleTree<Fq> = MerkleTree::new(depth);
        let new_data = Fq::from(10);

        merkle_tree.update_leaf(0, new_data, false).unwrap();
        merkle_tree.update_leaf(0, Fq::from(0), true).unwrap();

        assert_eq!(merkle_tree.leaves[0], Fq::from(0));

        let zero_hash = Fq::from(0);
        let hash_1 = MerkleTree::<Fq>::hash_pair(zero_hash, zero_hash);
        let hash_2 = MerkleTree::<Fq>::hash_pair(hash_1, hash_1);

        assert_eq!(merkle_tree.get_root_hash(), hash_2);
    }

    #[test]
    fn test_proof_and_verify() {
        let depth = 3;
        let mut merkle_tree: MerkleTree<Fq> = MerkleTree::new(depth);
        let new_data = Fq::from(10);

        merkle_tree.update_leaf(0, new_data, false).unwrap();

        let proof = merkle_tree.create_proof(new_data, 0).unwrap();
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
            proof: vec![fake_proof_data; depth],
        };

        let is_verified = merkle_tree.verify(invalid_proof);

        assert_eq!(is_verified, false);
    }

    #[test]
    fn test_create_proof_invalid_data() {
        let depth = 2;
        let mut merkle_tree: MerkleTree<Fq> = MerkleTree::new(depth);
        let new_data = Fq::from(10);

        merkle_tree.update_leaf(0, new_data, false).unwrap();

        let wrong_data = Fq::from(20);
        let result = merkle_tree.create_proof(wrong_data, 0);

        assert!(result.is_err());
    }

    #[test]
    fn test_new_with_inputs() {
        let depth = 2;
        let inputs = vec![Fq::from(1), Fq::from(2), Fq::from(3)];

        let merkle_tree = MerkleTree::<Fq>::new_with_inputs(depth, inputs.clone()).unwrap();

        assert_eq!(merkle_tree.leaves.len(), 4);
        assert_eq!(merkle_tree.tree.len(), depth);
        assert_eq!(merkle_tree.tree[0].len(), 2);
        assert_eq!(merkle_tree.tree[1].len(), 1);

        for (i, input) in inputs.iter().enumerate() {
            assert_eq!(
                merkle_tree.leaves[i],
                MerkleTree::<Fq>::compute_hash(*input)
            );
        }

        println!("tree is {:#?}", merkle_tree);

        assert_eq!(merkle_tree.leaves[3], Fq::from(0));

        let too_many_inputs = vec![Fq::from(1); 5];
        assert!(MerkleTree::<Fq>::new_with_inputs(depth, too_many_inputs).is_err());
    }
}
