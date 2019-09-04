# Data analysis, model selection and model calibration for "Hidden long-range memories..." by Kuchen, Becker et al.

This repository collects software for statistical analysis and model selection
used in the article "Hidden long-range memories of growth and cycle speed
correlate cell cycles in lineage trees", by E. Kuchen, N. Becker, N. Claudino
and T. Höfer. The software was written by E. Kuchen and N. Becker.

The repository is divided in three directories, relating to the Bifurcating
AutoRegressive Hidden Markov Model (bar-hmm/), to analysis of correlations in
presence of censoring and potential confounding variables (data-analysis/) and
to fitting of the mechanistic growth-progression model (growth-progression/),
respectively.

## Bifurcating AutoRegressive Hidden Markov Model (bar-hmm/)

### Installation and prerequisites

The code is written mainly in [OCaml](www.ocaml.org). It has been tested with
OCaml 4.07.1 and should run on both OS X and Linux. OCaml and the required
libraries are best installed via the [opam package manager](opam.ocaml.org).
The following OCaml libraries are required (listed with versions):

package        | tested version | prerequisite
:---           | ---:           | ---:
containers     | 2.6            |
gen            | 0.5.2          |
sequence       | 1.1            |
ctypes         | 0.15.1         |
ctypes-foreign | 0.4.0          |
csv            | 2.2            |
gsl            | 1.24.0         | GSL
lacaml         | 11.0.4         | OpenBLAS
npy            | 0.0.8          |
parmap         | 1.0-rc10       |
ocamlbuild     | 0.14.0         |

On OS X, the prerequisites can be installed via the homebrew package manager.
On Linux they should be available via the system package manager.

The tree data structure is based on the RoseTree module by S. Cruanes, licensed
as described in the file RoseTree.mli

The multivariate normal distribution integration routines in Fortran, adapted
from Alan Genz and Frank Bretz are included in the subdirectory
multivariate_normal_int/ .  These have to be compiled with a Fortran compiler
as a shared library and installed, as described in that directory, before
building the main program.

The tool ocamlbuild is used for building the main program. From the source
directory, to build the native code executable for evidence calculation, type
"ocamlbuild model_evidence.native". To run the program, type
"./model_evidence.native". To clean up the installation, type "ocamlbuild
-clean".

### Scope

The central piece of code in this directory is a library (hmtree.ml) for
evaluating the exact likelihood for a bifurcating autoregressive cell-cycle
model with a general number d of latent Gaussian variables. The model is
described in the article, and a prototype is included in the form of a
Mathematica notebook (recursive-likelihood-new.nb).

The library is used for evaluating model evidences for various
specialized versions of the BAR model, where different numbers of and couplings
between hidden variables are allowed (in model_evidence.ml).

The remaining files are supporting libraries and a submit script for PBS
computing clusters.


### List of relevant top-level files

- multivariate_normal_int/
	
	Multinormal integration

- tests/

	Basic internal tests

- tree_data/

	Raw data from experiments in .csv format. Copy to $HOME/tmp/ before use.

- export_tree_mma.ml

	Conversion of trees to Mathematica format for visualization

- hmtree.ml

	Main module for BAR-HMM likelihood evaluation

- import_csv.ml

	Import raw data from csv

- integrate_mvnd.ml

	Wrapper for multinormal integration

- models.ml

	Definition of variants of the BAR model

- model_evidence.ml

	Main driver script to calculate model evidences

- model_relcorr_evidence.ml

	Approximate model evidences based on correlations only

- recursive-likelihood-new.nb

	Prototype of the recursice BAR-HMM likelihood

- roseTree.ml

	Tree handling library

- roseTree.mli

	Interface thereof

- submit_evidence.ml

	Script to submit evidence calculations on a PBS batch server

- test_lib.ml

	Cell-cycle utility library

- tree.ml

	Tree-related utility library

- util.ml

	Basic utility library

- _tags

	Build definition file used by ocamlbuild



## Data analysis (data-analysis/)


Cell lineage tree data analysis pipeline.

The lineage tree data was analysed using MATLAB (version R2016b) and its
Statistics and Machine Learning toolbox. Some statistical tests were performed
in R (version 3.4.3) using the stats package.

To re-analyse the data, load MS_KuchenBeckeretal_datanalysis.m. This script
loads in the lineage tree data (published as supplemental data with the
article). The code is divided into code blocks, which can be run consecutively. Each block generates one figure panel presented in the manuscript. The code blocks are titled by the figure panels they generate.  

First the lineage tree control data (published as supplemental data: LineageTrees_controlReplicates_Neuroblastoma_KuchenBeckeretal) is loaded and a bias correction for finite observation time is performed. The corrected lineage trees are saved as ‘.mat’ files in the working directory. The corrected lineage trees are then again loaded and cycle length correlations and distributions analysed, within the whole population and with respect to the cell lineage, time and space.  


### List of files

- MS_KuchenBeckeretal_datanalysis.m

	Main data analysis file that calls the functions below

- PlotCorrelations.m

	Graphical presentation of correlation patterns

- TemporalDrift_PartialCorrelations.m

	Test for spatial proximity and temporal drift as potential causes for
correlations

- calcCorrSpearman.m

	Support file for correlation coefficient calculation

- dobootstrapTrees.m

	Bootstrapping for confidence intervals


## Growth-progression model (growth-progression/)

The growth-progression model is implemented in R (version 3.4.3) and simulates
cell cycle lengths within lineage trees as described in Kuchen, Becker et al.
The implementation does not require additional R packages to run. However, if output files should be saved for downstream analysis with MATLAB, package ‘R.matlab’ is required.  

Model, simulation and output parameters can be specified in the initial code section marked USER INPUT. Current settings present the settings used to simulate the results shown for rep2. Prototypic output figures are generated. The more refined figure panel versions shown in the publication were generated with MATLAB.

The user can decide between two different cell growth models, either exponential (used for rep, rap and esc cells) or logistic growth (used for –myc cells). 

The implemented model can simulate either full branched-out lineage trees (using checkstationary=0) or with checkstationary=1 full branched-out lineage trees up to a fixed number of generations (determined by ‘xgen’) with cells being randomly removed and lineage trees simulated for another ‘genssimu-xgen’ generations thereafter to avoid memory issues do to an excessively large data frame. The former (checkstationary=0) was used for the simulation results shown. Checkstationary=1 was only used to investigate the stationarity in cell age and size over time. 

The growth-progression model simulation is additionally provided as a
functionally equivalent but faster OCaml program. Both R and Ocaml versions
were tested for equivalence and used in various stages of the project.

Finally, lineage tree simulations were further analysed using MATLAB, file MS_KuchenBeckeretal_modelanalysis.m. This script performs a best fit parameter estimation using Approximate Bayesian Computing using simulated lineage trees under many different parameter combinations as input. In addition, the distribution of cell age and cell size over time can be analysed using simulation runs under checkstationary=1 (see above) as input. 


### List of files

- GrowthProgressionModel.R

	Model, simulation and output parameters can be specified in the initial code
section marked USER INPUT. Current settings present the settings used to
simulate the results shown for rep2. Prototypic output figures are generated.
(The more refined figure panel versions in the article were generated with
MATLAB.)


- ABC_parameterBounds_weighted_alldatasets.m

	Support file for ABC

- MS_KuchenBeckeretal_modelanalysis.m

	ABC parameter estimation for the growth-progression model, based on the set
of cycle correlations between related cells. Investigation of cell age and size homeostasis.

- simulate_growth_progression.ml

	Sampling (stochastic simulation) of the growth-progression model, Ocaml
version. Required libraries and building as described in the BAR model section.
To build, this file should be copied to the BAR model directory since it relies
on supporting libraries therein.

- submit_sim_gp.ml

	Script (standalone) to submit growth-progression simulations on a PBS batch
server.
