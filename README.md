# Bifurcating autoregressive hidden Markov model for cell cycle inheritance

This repository collects software for statistical analysis and model selection
"Hidden long-range memories of growth and cycle speed correlate cell cycles in
lineage trees", by E. Kuchen, N. Becker et al.

## Installation and prerequisites

The code is written mainly in [OCaml](www.ocaml.org). It has been tested with
OCaml 4.07.1 and should run on both OS X and Linux. OCaml and required
libraries are best installed via the [opam package manager](opam.ocaml.org).
The following OCaml libraries are required (listed with tested versions):

name    containers
version 2.6
name    csv
version 2.2
name    ctypes
version 0.15.1
name    ctypes-foreign
version 0.4.0
name    gen
version 0.5.2
name    gsl
version 1.24.0
name    lacaml
version 11.0.4
name    npy
version 0.0.8
name    parmap
version 1.0-rc10
name    sequence
version 1.1

The libraries lacaml and gsl reqiure a system-wide installation of openblas and
of the GSL. On OS X these can be installed via the homebrew package manager.

For building, the package ocamlbuild is needed in addition. From the source
directory, to build the native code executable for evidence calculation, type
"ocamlbuild model_evidence.native"
To clean the installation, type
"ocamlbuild -clean"

The tree handling is based on the RoseTree module by S. Cruanes, licensed as
described in the file RoseTree.mli

The multivariate normal distribution integration routines in Fortran, adapted
from versions by Alan Genz and Frank Bretz are included in the subdirectory
multivariate_normal_int/ .  These have to be compiled with a Fortran compiler
as a shared library and installed, as described in that directory, before
building the main program.


## Scope

This repository contains a library for evaluating the exact likelihood for a
bifurcating autoregressive cell-cycle model with d latent Gaussian variables.
The library is used for evaluating model evidences for various specialized
versions of the BAR model, where different numbers of and couplings between
hidden variables are allowed. This model is described in the article and
prototyped in a Mathematica notebook contained herein.

In addition, there are numerical simulations for generating trees, based on
both the BAR model and the growth-progression model as described in the paper.


## List of relevant top-level files

- multivariate_normal_int/
	
	Multinormal integration

- tests/

	basic internal tests

- tree_data/

	raw data from experiments in .csv format. copy to $HOME/tmp/ before use.

- export_tree_mma.ml

	conversion of trees to Mathematica format for visualization

- hmtree.ml

	main module for BAR-HMM likelihood evaluation	

- import_csv.ml

	import raw data from csv

- integrate_mvnd.ml

	wrapper for multinormal integration

- models.ml

	definition of variants of the BAR model

- model_evidence.ml

	main driver script to calculate model evidences

- model_relcorr_evidence.ml

	approximate model evidences based on correlations only

- recursive-likelihood-new.nb

	prototype of the recursice BAR-HMM likelihood

- roseTree.ml

	tree handling library

- roseTree.mli

	interface thereof

- simulate_growth_progression.ml

	sampling (stochastic simulation) of the growth-progression model	

- submit_evidence.ml

	script to submit evidence calculations on a PBS batch server

- submit_sim_gp.ml

	script to submit growth-progression simulations on a PBS batch server

- test_lib.ml

	cell-cycle utility library

- tree.ml

	tree-related utility library

- util.ml

	basic utility library

- _tags

	build definition file used by ocamlbuild
