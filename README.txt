This code repository contains the files required to reproduce the analysis and figures in the paper "Capacity Limits Lead to Information Bottlenecks for Ongoing Rapid Motor Behaviour."

Unless otherwise noted, all code is original to the analysis for these works.

1. DATA ANALYSIS
- The script automateRunningData.m goes through all four task types and calls the function processMultipleOHFiles.m for each block of trials. At the end, the results from each task's blocks are combined into a single .mat file for further analysis.

- The script characterizeAllTrialsByTaskType.m goes through all four task types and analyses subjects' performance in all task trials. A high-level, individual task description is provided.

- The script analyseSubjectPerformanceAcrossTrialTypes.m goes through four comparison: OH-OHA, TOH-TOHA, OH-TOH, and TOH-OHA. For each comparison, the performance of subjects who have performed both tasks is compared. This is generally done by looking at paired trials between the tasks.

NOTE: The KINARM utility functions in this project are freely available, with registration, from https://kinarm.com/support/software-downloads/

analyseSubjectPerformanceAcrossTrialTypes.m
automateRunningData.m
characterizeAllTrialsByTaskType.m
calculateScreenBounds.m
characterizeSingleTaskType.m
copyMetaData.m
getBasicStatistics.m
getBinAccuracy.m
getDecisionList.m
getEventsList.m
getStatesList.m
getSteadyStateInfo.m
getXYbyChannel.m
KINARM_add_hand_kinematics.m *- from KINARM
processMultipleOHFiles.m
processSingleOHFile.m
sort_trials.m *- from KINARM
zip_load.m *- from KINARM
ohObject.m
screenBounds.m
trialObject.m


2. BAYESIAN STATISTICAL ANALYSIS
- The script ssrDifferences.R performs Bayesian data analysis to parameterize the distribution of differences between steady-state rates for two different tasks.

- The script ssrLinearRegression.R performs Bayesian linear regression to parameterize the distribution of steady-state rates in one task using a linear function of another task's steady-state rates.

NOTE: These scripts rely on code written to accompany the book "Kruschke, J. K. (2015). Doing Bayesian Data Analysis, Second Edition: A Tutorial with R, JAGS, and Stan. Academic Press / Elsevier." These files (and many others) can be found at the book's website: https://sites.google.com/site/doingbayesiandataanalysis/software-installation.

DBDA2E-utilities.R
Jags-Ymet-Xmet-Mrobust-OHOHA.R
JAGS-Ymet-Xnom1gp-Mnormal-OHOHA.R
ssrDifferences.R
ssrLinearRegression.R


3. SSS SIMULATION
- The script SSSExperimentFramework.m parameterizes and runs an experiment to compare the performance of three kinds of simple serial systems (e.g., limited-capacity, parallel-serial processing systems). Simulation results are stored in a .mat file for further analysis.

SSSExperimentFramework.m
AdditiveWorkloadSerialSystem.m
FlexibleCapacitySerialSystem.m
MultiplicativeWorkloadSerialSystem.m


4. FIGURE PRODUCTION
- The script informationBottleneckPaperFigures.m takes the .mat files produced in DATA ANALYSIS and SSS SIMULATION and generates the paper's panels. In some cases, additional panels are produced, for example when only a single task is illustrated in the paper it is generally possible to generate the same panel for the other three tasks as well.

informationBottleneckPaperFigures.m


5. UTILITIES
- A collection of functions that are called throughout the other scripts.

alignTrials.m
blandAltmanPlot.m
cividis.m *- from MATLAB File Exchange (https://www.mathworks.com/matlabcentral/fileexchange/66493-alexhenderson-cividis)
dir2.m
matFileCombiner.m
mergeStructs.m
savefigas.m
scatterAndRegressionWithCIs.m
