# Information-Bottlenecks-in-Rapid-Motor-Behaviour
Analysis and figure generation code for the paper "Capacity Limits Lead to Information Bottlenecks for Ongoing Rapid Motor Behaviour."

## Note about licenses

The majority of this paper's code is available under the GPL-3.0 license. There are two broad exceptions to this.

Code derived from John Kruschke's book "Doing Bayesian Data Analysis, Second Edition: A Tutorial with R, JAGS, and Stan" is freely available from the book's [website](https://sites.google.com/site/doingbayesiandataanalysis/software-installation).

Code from Kinarm is copyright BKIN Technologies Ltd and is shared here with permission. These files are included for the purposes of reproducibility. The latest versions are freely available, with registration, from the [Kinarm website](https://kinarm.com/support/software-downloads/).

## Code

1. Data Analysis

Scripts and functions for analysing raw .c3d files produced by a Kinarm Exoskeleton or Endpoint robot. Calculates a number of trial-level statistics for the OH, OHA, TOH, and TOHA tasks.

NOTE: The Kinarm utility functions in this project are freely available, with registration, from the [Kinarm website](https://kinarm.com/support/software-downloads/).

2. Bayesian Statistical Analysis

Scripts for performing Bayesian statistical analysis.

NOTE: These scripts rely on code written to accompany the book "Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier." These files (and many others) can be found at the book's [website](https://sites.google.com/site/doingbayesiandataanalysis/software-installation).

3. SSS Simulation

Parameterize and run experiments to compare the performance of three kinds of simple serial systems (e.g., limited-capacity, parallel-serial processing systems).

4. Figure Production

Generates panels for the paper's figures. In some cases, additional panels are produced. For example, when only a single task is illustrated in the paper it is generally possible to generate the same panel for the other three tasks as well.

5. Utilities

A collection of functions that are called throughout the other scripts.

## Documentation

1. README.txt

An in-depth listing and description of the scripts and functions contained in each folder.
