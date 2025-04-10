# ZK Research Implementations

A Rust implementation of various zero-knowledge proof protocols and cryptographic primitives.

## Overview

This library provides implementations of several key cryptographic components for zero-knowledge proofs:

- GKR Protocol (Goldwasser-Kalai-Rothblum)
- KZG Polynomial Commitment Scheme
- Fiat-Shamir Transform
- Shamir Secret Sharing
- Multilinear and Univariate Polynomial Operations
- Fast Fourier Transform (FFT)
- Merkle Tree
- Sum-Check Protocol

## Components

### GKR Protocol

The GKR protocol implementation consists of:

- **Circuit Representation** (`gkr_circuit.rs`):

  - Defines arithmetic circuits with add and multiply gates
  - Supports layer-wise circuit evaluation
  - Represents computation as multilinear polynomials

- **Protocol Implementation** (`gkr_protocol.rs`):
  - Provides proof generation and verification
  - Implements the sum-check protocol for GKR
  - Handles polynomial operations for efficient verification

### Polynomial Commitment Schemes

- **KZG Commitment** (`kzg.rs`):
  - Kate-Zaverucha-Goldberg polynomial commitment scheme
  - Supports trusted setup generation
  - Implements commitment, opening and verification
  - Uses BLS12-381 pairing-friendly curve

### Fiat-Shamir Transform

- **Transcript Generation** (`fiat_shamir_transcript.rs`):
  - Implements non-interactive proof generation
  - Provides secure random challenge derivation
  - Uses Keccak256 for hashing

### Shamir Secret Sharing

- **Secret Sharing** (`shamir_secret_sharing.rs`):
  - Implements threshold-based secret sharing
  - Supports secret recovery from shares
  - Uses polynomial interpolation techniques

### Polynomial Operations

- **Multilinear Polynomials** (`multilinear_polynomial_evaluation.rs`):

  - Implements multilinear polynomial representation and operations
  - Supports efficient evaluation and partial evaluation
  - Provides basic arithmetic operations (add, multiply, subtract)

- **Composite Polynomials** (`composed_polynomial.rs`):

  - Implements product and sum polynomial structures
  - Supports operations on complex polynomial compositions

- **Univariate Polynomials** (`univariate_polynomial_dense.rs`):
  - Dense representation of univariate polynomials
  - Implements polynomial interpolation
  - Provides evaluation and arithmetic operations

### Fast Fourier Transform (FFT)

- **FFT Implementation** (`fft.rs`):
  - Implements Discrete Fourier Transform (DFT)
  - Provides polynomial evaluation via FFT
  - Supports polynomial interpolation via inverse FFT
  - Optimized for finite field operations

### Merkle Tree

- **Merkle Tree** (`merkle_tree.rs`):
  - Implements a binary Merkle tree structure
  - Supports proof generation and verification
  - Provides leaf updates and path recomputation
  - Uses field elements as leaf values

### Sum-Check Protocol

- **Sum-Check Protocol** (`sum_check_protocol.rs`):
  - Implements the sum-check protocol for multilinear polynomials
  - Provides proof generation and verification
  - Includes specialized GKR protocol integration
  - Supports composed polynomial structures

### Sample Tests

- **Fibonacci Evaluation** (`fibonacci_evaluation.rs`):
  - Demonstrates polynomial interpolation with Fibonacci sequence
  - Shows practical application of univariate polynomials
  - Provides verification of Fibonacci properties

## Technical Details

### GKR Circuit

The GKR circuit implementation allows for creating and evaluating arithmetic circuits with:

- Support for addition and multiplication gates
- Layer-based circuit structure
- Automatic conversion to multilinear polynomial representation

### KZG Polynomial Commitments

The KZG scheme enables:

- Creating commitments to polynomials
- Opening commitments at specific points
- Verifying openings using bilinear pairings
- Batched operations for efficiency

### Fiat-Shamir Transcript

The transcript system:

- Ensures non-interactive zero-knowledge proofs
- Provides deterministic challenge generation
- Maintains transcript state throughout proof generation

### Multilinear Polynomials

The multilinear polynomial implementation provides:

- Evaluation on binary cube vertices
- Efficient partial evaluation algorithms
- Operations compatible with the GKR protocol

### Fast Fourier Transform

The FFT implementation provides:

- Efficient polynomial evaluation in O(n log n) time
- Polynomial interpolation via inverse FFT
- Optimized for finite field operations
- Support for various field types

### Merkle Tree

The Merkle tree implementation offers:

- Efficient membership proofs
- Logarithmic verification complexity
- Support for dynamic updates
- Integration with field elements

### Sum-Check Protocol

The sum-check protocol implementation provides:

- Interactive proof system for sum verification
- Non-interactive variant via Fiat-Shamir
- Integration with GKR protocol
- Support for composed polynomial structures

## Dependencies

- `ark-ff`: Finite field operations
- `ark-ec`: Elliptic curve operations
- `ark-bls12-381`: BLS12-381 curve implementation
- `ark-bn254`: BN254 curve implementation
- `ark-poly`: Polynomial operations
- `sha3`: Keccak256 hashing
- `rand`: Random number generation

## Usage

All components are made to be easy to understand and use. The library provides a comprehensive set of tools for implementing zero-knowledge proofs and cryptographic primitives.

## Testing

Each component includes comprehensive unit tests demonstrating functionality and correctness. Run tests using:

```bash
cargo test
```

## License

[License information would go here]
