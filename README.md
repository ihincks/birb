# Code for "Bayesian Inference for Randomized Benchmarking Protocols"

- Ian Hincks
- Joel J. Wallman
- Chris Ferrie
- Chris Granade
- David G. Cory

## Introduction

This repository contains all source code necessary to reproduce the results found in the paper _Bayesian Inference for Randomized Benchmarking Protocols_. Plase submit questions about using this repository as GitHub issues.

## Requirements

Both Mathematica and Python are used.

#### Mathematica

File extensions `.nb` are to be run with _Wolfram Mathematica 11_ (though earlier and later versions will likely work too). As this is paid software, we have also made the notebooks available in the `.cdf` format where they can be viewed, but not run, using the free [Wolfram CDF Player](https://www.wolfram.com/cdf-player/). 

The following third party _Mathematica_ packages are required:

 - [QuantumUtils for Mathematica](https://github.com/QuantumUtils/quantum-utils-mathematica)
 - [MaTeX](https://github.com/szhorvat/MaTeX)

#### Python

We recommend using a conda virtual environment, installed with 

```bash
$ conda install nb_conda # Optional, exposes new environment to Jupyter.
$ conda env create -f environment.yml
```

These requirements can be installed automatically with pip:

```bash
$ pip install -r requirements.txt
```

Or with pip and virtualenv:
```
$ virtualenv env/
$ env/scripts/activate.sh # Use ".ps1" instead of ".sh" on Windows.
$ pip install -r requirements.txt
```
## Index

#### By Figure Number

- [Figure 1](https://github.com/ihincks/birb/blob/master/fig/betabin-crb.pdf): 
`src/sequence-reuse-calculations.*`, subsection _"Weighted Cramer-Rao Bound Plots"_
- [Figure 2](https://github.com/ihincks/birb/blob/master/fig/unitarity-optimal-switching-cost.pdf): 
`src/sequence-reuse-calculations.*`, subsection _"Generalization to Finite Switching Cost"_
- [Figure 3](https://github.com/ihincks/birb/blob/master/fig/different-noise-models-summary.pdf): 
` `, 
- [Figure 4](https://github.com/ihincks/birb/blob/master/fig/different-noise-models-comparison.pdf): 
` `, 
- [Figure 5](https://github.com/ihincks/birb/blob/master/fig/different-noise-models-survival-dists.pdf): 
` `, 
- [Figure 6](https://github.com/ihincks/birb/blob/master/fig/low-data.pdf): 
` `, 
- [Figure 7](https://github.com/ihincks/birb/blob/master/fig/pathological-survival-distributions.pdf): 
` `, 
- [Figure 8](https://github.com/ihincks/birb/blob/master/fig/lrb-posterior.pdf): 
` `, 

#### By Table Number

- Table 2: `src/beta-reparameterizations.*`

#### By Equation Number
