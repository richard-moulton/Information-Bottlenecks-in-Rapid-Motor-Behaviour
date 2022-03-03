# Performing Bayesian Linear Regression for Steady-State Rates from the OH family of tasks
# Based on code and analysis from "Doing Bayesian Data Analysis, 2nd Edition",
# John Kruschke, 2014.
# Copyright (C) 2022  Richard Hugh Moulton
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# 3 March, 2022

# DO NOT EXECUTE CODE ABOVE HERE. R's MEMORY IS CLEARED.

# Optional preliminaries.
graphics.off() # Closes all of R's graphics windows.
rm(list=ls())  # Clears all of R's memory.

# Ensure libraries and functions are appropriately loaded
library(ggplot2)
theme_set(theme_classic())
library(gridExtra)
library(stats)
library(plyr)
library(MASS)
source("assessNormality.R")

# Load the the necessary code from Bayesian Data Analysis
source("DBDA2E-utilities.R")
source("Jags-Ymet-Xmet-Mrobust-OHOHA.R")

# Set global variables for figures
graphFileType = "png"

# Determine which comparison to make
# Comp Flag
# 1: OH-OHA, 2: OH-TOH, 3: TOH-TOHA, 4: TOH-OHA

for (compFlag in 4:4)
{
  # Load the data file 
  fileToLoad = switch(compFlag,"SubjectSteadyStateRates(OH-OHA-SingleTrial).csv",
                      "SubjectSteadyStateRates(OH-TOH-SingleTrial).csv",
                      "SubjectSteadyStateRates(TOH-TOHA-SingleTrial).csv",
                      "SubjectSteadyStateRates(TOH-OHA-SingleTrial).csv");
  dataFrame = read.csv( file=fileToLoad )
  trialShortNames = dimnames(dataFrame)[[2]]
  trialLongNames = c(switch(trialShortNames[1],"OH"="Object-Hit","OHA"="Object-Hit-and-Avoid",
                            "TOH"="Turbo Object-Hit","TOHA"="Turbo Object-Hit-and-Avoid"),
                     switch(trialShortNames[2],"OH"="Object-Hit","OHA"="Object-Hit-and-Avoid",
                            "TOH"="Turbo Object-Hit","TOHA"="Turbo Object-Hit-and-Avoid"));
  ssr1 = dataFrame[,1]
  ssr2 = dataFrame[,2]
  fileNameRoot = switch(compFlag,"InformationBottleneck-SSRRegression-OHOHA-",
                        "InformationBottleneck-SSRRegression-OHTOH-",
                        "InformationBottleneck-SSRRegression-TOHTOHA-",
                        "InformationBottleneck-SSRRegression-TOHOHA-");
 
  #------------------------------------------------------------------------------- 
  # Generate the MCMC chain:
  mcmcCoda = genMCMC( data=dataFrame , xName=trialShortNames[1] , yName=trialShortNames[2] , 
                      numSavedSteps=20000 , saveName=fileNameRoot )
  #------------------------------------------------------------------------------- 
  # Display diagnostics of chain, for specified parameters:
  parameterNames = varnames(mcmcCoda) # get all parameter names
  for ( parName in parameterNames ) {
    diagMCMC( codaObject=mcmcCoda , parName=parName , 
              saveName=fileNameRoot , saveType=graphFileType )
  }
  # Get summary statistics of chain:
  summaryInfo = smryMCMC( mcmcCoda , 
                          compValBeta0=0.0, ropeBeta0=c(-0.12, 0.12),
                          compValBeta1=1.0 , ropeBeta1=c(0.95,1.05) ,
                          saveName=fileNameRoot )
  # Display posterior information:
  plotMCMC( mcmcCoda , data=dataFrame , xName=trialShortNames[1] , yName=trialShortNames[2] , 
            pairsPlot=TRUE , showCurve=FALSE ,
            saveName=fileNameRoot , saveType=graphFileType )
  #------------------------------------------------------------------------------- 
  
}