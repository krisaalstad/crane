<h1 align="center">crane</h1>
<p align="center">
  <em>Cryospheric Retrospective Analysis using Nested Ensembles</em>
</p>

<p align="center">
  <img src=".github/assets/logo.png" alt="CRANE Logo" width="220" />
</p>

<p align="center">
  <em>crane is built on the back of a <a href="https://en.wikipedia.org/wiki/Turtles_all_the_way_down">two-turtle</a> hierarchical Bayesian model.</em>
</p>


**crane** is the code repository accompanying the manuscript *Climate-informed cryospheric reanalysis via hierarchical Bayesian data assimilation*.

All the source code is in the subdirectory `src`, while the run scripts are `mainsnow.m` and `mainglacier.m` for the snow and glacier experiments, respectively. These scripts set up some of their own custom run parameters, while the rest are defined in `setpars.m`. The main scripts are configured to run as is, but users are free to change settings by toggling switches. For example, setting `p.assimF=1` and `p.assimD=0` in the snow script switches from the SWE data assimilation experiment to the FSCA data assimilation experiment.

The source code in `src` consists of several interlinked functions. The key functions are: the **cryospheric reanalysis** routine (`CRA.m`); nearly identical degree-day (temperature-index) **models** for glaciers (`DDMg.m`) and snow (`DDMs.m`); an **ensemble Kalman** analysis (`EnKA.m`); the Multiple Adaptively Guided Particle-adjusted Iterative Ensemble Smoother (**MAGPIES**, `MAGPIES.m`); an **Expectation Maximization** routine (`particleEM.m`); **Particle Markov Chain Monte Carlo** using the Robust Adaptive Metropolis algorithm (`RAMP.m`); a routine to recycle MAGPIES for new hyperparameters (`reMAGPIES.m`); and, last but not least, an **Adaptive Particle Batch Smoother** routine for hyperparameter inference using the weight-based evidence approximation for Adaptive Multiple Importance Sampling (`WAMIS.m`). All of these key functions are described in the corresponding sections of the manuscript.

Running the code requires a MATLAB programming environment (tested with R2021a). For convenience, all the input data required to run the experiments are downloaded automatically the first time one of the main scripts is run. Without Markov chains, running an entire experiment for all sites takes only a few minutes on a standard laptop. To run Markov chains (`p.doMCMC=1`) as part of the PMCMC benchmarking for all experiments on wall-clock timescales of days rather than weeks or months, you will also need access to 10 or more cores, ideally on an external server. To avoid having to rerun any experiments, including the Markov chains, all the output data analyzed in the manuscript are downloaded automatically the first time one of the plotting scripts (`pardynfig.m`, `plot_evaluate.m`, and `timefigs.m`) or the evaluation script (`evaluate.m`) is run. The input and output data are fetched from the accompanying [Zenodo](https://doi.org/10.5281/zenodo.22938686) dataset.

Bug reports are welcome as [GitHub issues](https://github.com/krisaalstad/crane/issues). For other questions, or if you would like help to adapt **crane** for your own work, feel free to get in touch by email: kristoffer.aalstad@geo.uio.no

*Disclaimer*: The code in **crane** is implemented entirely in MATLAB using elementary array programming and vectorization, which makes it fast. In principle, it should be straightforward to translate to other languages such as Python or Julia, especially with the aid of LLM-based coding agents. The code is entirely human-written, but parts of it have been screened for bugs by LLMs. The corresponding author takes full responsibility for any mistakes in the code and would be happy to be informed of any bugs. A disadvantage of MATLAB is that it is proprietary commercial software, which many researchers both inside and outside academia do not have access to. An advantage of MATLAB is that it is proprietary commercial software, which means a well-maintained environment, self-contained code, highly optimized array operations, and few low-level bugs. While MATLAB is not open source, it is source-available, so much of its lower-level code, such as native `*.m` functions, can still be inspected by users. We encourage other researchers to adapt and translate this code to other programming environments. For example, several of the routines in `src` are already available in the Python-based Multiple Snow Data Assimilation system [MuSA](https://github.com/ealonsogzl/MuSA).

*Logo*: The **crane** logo was generated with Nano Banana (Google Gemini) and edited by the author with GIMP.
