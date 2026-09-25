<h1 align="center">crane</h1>
<p align="center">
  <em>Cryospheric Retrospective Analysis using Nested Ensembles</em>
</p>

<p align="center">
  <img src=".github/assets/logo.png" alt="CRANE Logo" width="220" />
</p>

<p align="center">
  <em>Crane is built on the back of a <a href="https://en.wikipedia.org/wiki/Turtles_all_the_way_down">two-turtle</a> hierarchical Bayesian model.</em>
</p>




This is the code repository accompanying the manuscript *Climate-informed cryospheric reanalysis via hierarchical Bayesian data assimilation*. 

All the source code is in the subdirectory `src` while the run scripts are `mainsnow.m` and `mainglacier.m` for the snow and glacier experiments, respectively. These scripts setup some of their own custom parameters, while the remaining parameters are defined in `setpars.m`. These main scripts are configured to be run as is, but the user is free to switch settings by turning on and off switches, for example setting `p.assimF=1` and `p.assimD=0` in the snow script switches from the SWE data assimilation to the FSCA data assimilation experiment. 

The source code in `src` consists of several inter-linked functions, where the key functions are the: **cryospheric reanalysis** routine in `CRA.m`, nearly identical degree day (temperature) **models** for glaciers (`DDMg.m`) and snow (`DDMs.m`), an **ensemble Kalman** analysis (`EnKA.m`), Multiple Adaptively Guided Particle-adjusted Iterative Ensemble Smoother (**MAGPIES**, in `MAGPIES.m`), Expectation Maxmization routine (`particleEM.m`), **Particle Markov Chain Monte Carlo** using the Robust Adaptive Metropolis algorithm in (`RAMP.m`), a routine to recycle MAGPIES for new hyperparameters (`reMAGPIES.m`), and last but not least an **Adaptive Particle Batch Smoother** routine for hyperparaemter inference using the weight-based evidence approximation for Adaptive Multiple Importance Sampling (`WAMIS.m`). All of these key functions are described in the correspondingly **sections** in the manuscript.

Running the code requires access to a MATLAB (tested with R2021a) programming environment. For convenience, all the input data required to run the experiments are downloaded automatically the first time one of the main scripts is run. Without Markov chains, running an entire experiment for all sites takes a few minutes only on a standard laptop. To run Markov chains (`p.doMCMC=1`) as a part of PMCMC benchmarking for all experiments on wall-clock timescales of days and not weeks or months, you also need access to 10 (or more) cores ideally on an extenal server. To avoid having to rerun any experiments, including Markov chains, all the output data analyzed in the manuscript is automatically downloaded the first time one of the plotting scripts (`pardynfig.m`,`plot_evaluate.m`, and `timefigs.m`) or the evaluation script (`evaluate.m`) are run. The input and output data are fetched from the accompanying [Zenodo](https://doi.org/10.5281/zenodo.22938686) dataset. 

*Disclaimer*: The code in **crane** is implemented entirely in MATLAB, based on elementary array programming and vectorization, which makes it fast. In theory, it can also easily be translated to other programming languages such as Python or Julia, especially with the aid of LLM-based coding agents. The code is all human generated, but parts of it have also been screened for bugs by Large Language Models (LLMs). The corresponding author takes full responsibility for mistakes in the code, and would be happy to be informed about any bugs. A disadvantage with MATLAB is that it is proprietary commercial software, which is not available to many researchers both inside and outside of academia. An advantage with MATLAB is that it is proprietary software, so the code is self-contained, highly optimized for array programming, and unlikely to contain low level bugs. While MATLAB is not strictly open-source, it is source-available such that users can still inspect much of the lower level code, such as native `*.m` functions. We encourage other researchers to adapt and translate this code to other programming environments. For example, most of the routines in `src` are already available in the Python-based Multiple Snow Data Assimilation system [MuSA](https://github.com/ealonsogzl/MuSA).
