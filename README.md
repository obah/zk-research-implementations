# Zero Knowledge Proof Library

A Rust implementation of various zero-knowledge proof protocols and cryptographic primitives.

## Overview

This library provides implementations of several key cryptographic components for zero-knowledge proofs:

- GKR Protocol (Goldwasser-Kalai-Rothblum)
- KZG Polynomial Commitment Scheme
- Fiat-Shamir Transform
- Shamir Secret Sharing
- Multilinear and Univariate Polynomial Operations

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

## Dependencies

- `ark-ff`: Finite field operations
- `ark-ec`: Elliptic curve operations
- `ark-bls12-381`: BLS12-381 curve implementation
- `ark-bn254`: BN254 curve implementation
- `sha3`: Keccak256 hashing
- `rand`: Random number generation

## Usage

All components are made to be easy to understand and use.

## Testing

Each component includes comprehensive unit tests demonstrating functionality and correctness.
