# Federated ADMM-IRLS

This repository contains the implementation and mathematical derivation of a federated ADMM-IRLS framework for high-dimensional sparse logistic regression with binary outcomes across multiple data centers.

Each participating center retains its patient-level predictor matrix and binary outcome vector locally. Within each ADMM iteration, every center uses iteratively reweighted least squares (IRLS) to compute a local coefficient update from its own data. Only model-level quantities are communicated for global consensus estimation.

The local coefficient updates are coordinated through the Alternating Direction Method of Multipliers (ADMM). The central aggregation step applies soft-thresholding to obtain a common sparse coefficient vector, while the dual-variable updates enforce agreement across participating centers.

The framework enables collaborative sparse logistic-regression modeling, coefficient estimation, variable selection, and prediction without pooling or directly sharing patient-level data. Its privacy-preserving property is based on federated data locality rather than formal cryptographic privacy guarantees.

## Mathematical Derivation

The complete Federated ADMM-IRLS formulation, local IRLS derivation, global and dual updates, convergence diagnostics, and algorithm are available on the project webpage:

[View the Federated ADMM-IRLS](https://atikur616.github.io/Federated-ADMM-IRLS-for-Logistic-Regression/)

## Repository Files

- `index.html` — mathematical formulation, derivation, convergence diagnostics, and algorithm
- Federated ADMM-IRLS R code — implementation of the proposed method
- `README.md` — project overview
