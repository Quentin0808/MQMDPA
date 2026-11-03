# PhadSimu — MQMDPA MATLAB Implementation

## Overview

**PhadSimu** is a computational thermodynamics project developed by **Kun Wang**.

This repository provides a MATLAB implementation of the **Modified Quasichemical Model in the Distinguishable Pair Approximation (MQMDPA)** for calculating the thermodynamic properties of non-ideal solutions.

The MQMDPA describes manifold short-range ordering (MSRO) by distinguishing different types of nearest-neighbor pairs. It provides a thermodynamic framework for representing complex mixing behavior and local atomic interactions in multicomponent solutions.

The implementation is intended for research, model development, and numerical verification.

## Model Description

The MQMDPA accounts for manifold short-range ordering through distinguishable heteroatomic pairs and their associated Gibbs energy contributions.

The thermodynamic formulation accounts for the energies of all bond types, and the configurational entropy associated with their distribution.

Equilibrium pair populations are determined by minimizing the Gibbs energy subject to the relevant composition and pair-balance constraints.

## Numerical Implementation

The Modified Quasichemical Model in the Distinguishable Pair Approximation (MQMDPA) is implemented in MATLAB, with particular emphasis on the numerical treatment of manifold short-range ordering.

The equilibrium distribution of distinguishable pairs is determined through Gibbs energy minimization using the Broyden–Fletcher–Goldfarb–Shanno (BFGS) quasi-Newton algorithm.

The BFGS method approximates the Hessian matrix using gradient information, avoiding the need to explicitly calculate the full second-derivative matrix.

Key Features

Thermodynamic property calculations based on the MQMDPA.

Representation of multiple short-range ordering states through distinguishable pairs.

Gibbs energy minimization using the BFGS quasi-Newton method.

Numerical optimization routines independently implemented in MATLAB.

No dependency on MATLAB Optimization Toolbox functions.

The current implementation focuses on thermodynamic property calculations and the determination of equilibrium pair distributions. General multiphase equilibrium calculations are not included in this release.

## Requirements

The source code is written in MATLAB.

MATLAB may be required to execute the implementation. Compatibility with GNU Octave has not yet been verified.

The numerical algorithms are implemented without relying on MATLAB Optimization Toolbox routines.

## Usage

1. Clone or download this repository.
2. Open the MATLAB source files in a compatible environment.
3. Configure the model parameters and composition according to the supplied examples or input files.
4. Run the relevant calculation script.
5. Examine the calculated thermodynamic properties and equilibrium pair distributions.

Specific entry points, input formats, and examples should be identified in the source code documentation.

## Author

**Kun Wang**

Researcher in Computational Thermodynamics and Materials Modeling

Research interests include:

- Computational thermodynamics
- CALPHAD modeling
- Short-range ordering in solutions
- Thermodynamic model development
- Numerical methods for Gibbs energy minimization

## Citation and Acknowledgment

If you use this implementation in academic research, publications, or software development, please acknowledge the author's contribution and cite the corresponding software release or associated publications.

Suggested software citation:

Wang, K. *PhadSimu: MATLAB Implementation of the Modified Quasichemical Model in the Distinguishable Pair Approximation (MQMDPA).* GitHub repository.

https://github.com/Quentin0808/MQMDPA

Wang, K. et al., Extension of the modified quasichemical model in the distinguishable-pair approximation to multicomponent solutions via an on-the-fly interpolation framework, Calphad, 2025, 90, 102852. 

https://doi.org/10.1016/j.calphad.2025.102852


A version-specific DOI may be added in a future release.

## License

No open-source license has been specified for this repository at present.

The source code is publicly accessible, but public availability does not automatically grant permission for unrestricted reuse, modification, or redistribution.

Please contact the rights holder regarding permissions not otherwise granted.

## Disclaimer

This software is provided for research and educational evaluation. Users are responsible for independently verifying numerical results and determining suitability for their intended applications.