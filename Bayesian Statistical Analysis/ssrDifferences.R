# Bayesian Data Analysis on Steady State Rates in OH family of tasks
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

# Load the the necessary code from Bayesian Data Analysis
source("DBDA2E-utilities.R")
source("Jags-Ymet-Xnom1grp-Mnormal.R")

# Set global variables for figures
graphFileType = "png"

# Determine which comparison to make
# Comp Flag
# 1: OH-OHA, 2: OH-TOH, 3: TOH-TOHA, 4: TOH-OHA
# Population Flag
# 1: full population
# 0: quartile breakdown
populationFlag = FALSE;

for (compFlag in 1:4)
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
  if (populationFlag)
  {
    quartileFlag = 'P'
    ssr1 = dataFrame[,1]
    ssr2 = dataFrame[,2]
    fileNameRoot = switch(compFlag,"InformationBottleneck-SteadyStateRates-OHOHA-",
                          "InformationBottleneck-SteadyStateRates-OHTOH-",
                          "InformationBottleneck-SteadyStateRates-TOHTOHA-",
                          "InformationBottleneck-SteadyStateRates-TOHOHA-");
    deltaSSR = ssr2 - ssr1
    
    # Determine appropriate ROPEs to use for each parameter
    muROPE = switch(trialShortNames[1],"OH"=c(-0.3,0.3),"TOH"=c(-0.48,0.48))
    sigmaCOMP = sqrt(sd(ssr1)^2 + sd(ssr2)^2)
    sigmaROPE = c(0.9*sigmaCOMP,1.1*sigmaCOMP)
    esROPE = c(-1,1)
    
    # How precise can we be about the parameter theta, the reduction in steady-state
    # rate that occurs when a subject moves from one task to another?
    # From Jags-Ymet-Xnom1grp-Mnormal-Example.R 
    #------------------------------------------------------------------------------- 
    #------------------------------------------------------------------------------- 
    # Generate the MCMC chain:
    mcmcCoda = genMCMC( data=deltaSSR , numSavedSteps=20000 , saveName=fileNameRoot )
    #------------------------------------------------------------------------------- 
    # Display diagnostics of chain, for specified parameters:
    parameterNames = varnames(mcmcCoda) # get all parameter names
    for ( parName in parameterNames ) {
      diagMCMC( codaObject=mcmcCoda , parName=parName , 
                saveName=fileNameRoot , saveType=graphFileType )
    }
    #------------------------------------------------------------------------------- 
    # Get summary statistics of chain:
    summaryInfo = smryMCMC( mcmcCoda , 
                            compValMu=0 , ropeMu=muROPE ,
                            compValSigma=1 , ropeSigma=sigmaROPE ,
                            compValEff=0.0 , ropeEff=esROPE ,
                            saveName=fileNameRoot )
    show(summaryInfo)
    write.csv(summaryInfo,file=paste0(fileNameRoot,"SummaryInfo-P.csv"))
    # Display posterior information:
    plotMCMC( mcmcCoda , data=deltaSSR , 
              compValMu=0 , ropeMu=muROPE ,
              compValSigma=sigmaCOMP , ropeSigma=sigmaROPE ,
              compValEff=0.0 , ropeEff=esROPE ,
              pairsPlot=TRUE , showCurve=FALSE ,
              saveName=fileNameRoot , saveType=graphFileType )
    #------------------------------------------------------------------------------- 
    
    # Create a graph overlaying the three probability density functions of interest.
    healthyControlMean = switch(trialShortNames[1],"OH"=0.0665,"TOH"=0.0665)
    healthyControlSD = switch(trialShortNames[1],"OH"=0.2716,"TOH"=0.10)
    
    set.seed(1234)
    numSamples = 100000
    df <- data.frame(
      pdf=factor(rep(c("1-Region of Practical Equivalence","2-Healthy Control Model","3-Posterior Distribution"), each=numSamples)),
      samples=c(runif(numSamples, min=muROPE[1], max=muROPE[2]),
                rnorm(numSamples, mean=healthyControlMean, sd=healthyControlSD),
                rnorm(numSamples, mean=summaryInfo[1,3], sd=summaryInfo[2,3]))
    )
    head(df)
    
    means <- ddply(df, "pdf", summarise, grp.mean=mean(samples))
    head(means)
    
    cols <- c("1-Region of Practical Equivalence" = "#58d68d", "2-Healthy Control Model" = "#495ab0", "3-Posterior Distribution" ="#e84133")
    
    ### GRAPH ONE
    # Use semi-transparent fill
    openGraph(height=10,width=15)
    p1<-ggplot(df[df$pdf %in% "1-Region of Practical Equivalence",], aes(x=samples, fill=pdf)) +
     coord_cartesian(xlim = c(-1.5, 1.5), ylim = c(0.0, 1.8)) +
     geom_area(aes(y = ..density..),stat ="bin", position="identity", alpha=0.4) +
     theme(legend.title = element_blank()) +
     labs(title = "Comparing Density Functions",
          subtitle = paste(trialLongNames[1]," vs. ",trialLongNames[2]),
          x = "Difference in Steady State Rates (Hz)", y = "Probability Density",
          colour = "Probability Density Function", fill = "Probability Density Function") + 
     scale_colour_manual(values = cols,aesthetics = c("colour", "fill"))
    p1
    
    # Add mean lines
    p1+geom_vline(data=means[means$pdf %in% "1-Region of Practical Equivalence",], aes(xintercept=grp.mean, color=pdf),
                 linetype="dashed")
    
    saveGraph(paste0(fileNameRoot,"PDFs-1"),graphFileType)
    
    ### GRAPH TWO
    # Use semi-transparent fill
    openGraph(height=10,width=15)
    p2<-ggplot(df[df$pdf %in% c("1-Region of Practical Equivalence","2-Healthy Control Model"),], aes(x=samples, fill=pdf)) +
     coord_cartesian(xlim = c(-1.5, 1.5), ylim = c(0.0, 1.8)) +
     geom_area(aes(y = ..density..),stat ="bin", position="identity", alpha=0.4) +
     theme(legend.title = element_blank()) +
     labs(title = "Comparing Density Functions",
          subtitle = paste(trialLongNames[1]," vs. ",trialLongNames[2]),
          x = "Difference in Steady State Rates (Hz)", y = "Probability Density",
          colour = "Probability Density Function", fill = "Probability Density Function") + 
     scale_colour_manual(values = cols,aesthetics = c("colour", "fill"))
    p2
    
    # Add mean lines
    p2+geom_vline(data=means[means$pdf %in% c("1-Region of Practical Equivalence","2-Healthy Control Model"),], aes(xintercept=grp.mean, color=pdf),
                linetype="dashed")
    
    saveGraph(paste0(fileNameRoot,"PDFs-2"),graphFileType)
    
    ### GRAPH THREE
    # Use semi-transparent fill
    openGraph(height=10,width=15)
    
    p3<-ggplot(df, aes(x=samples, fill=pdf)) +
      coord_cartesian(xlim = c(-1.5, 1.5), ylim = c(0.0, 1.8)) +
      geom_area(aes(y = ..density..),stat ="bin", position="identity", alpha=0.6) +
      theme(legend.position = "none", axis.text = element_text(size = 72), axis.line = element_line(size = 6)) +
      labs(x="",y="") + 
      scale_colour_manual(values = cols,aesthetics = c("colour", "fill"))
    
    show(p3)
    saveGraph(paste0(fileNameRoot,"PDFs-3-",quartileFlag),graphFileType)
    
  }
  else
  {
    for (quartileFlag in 1:4)
    {
      ssr1 = dataFrame[,1]
      ssr2 = dataFrame[,2]
      fileNameRoot = switch(compFlag,"InformationBottleneck-SteadyStateRates-OHOHA-",
                            "InformationBottleneck-SteadyStateRates-OHTOH-",
                            "InformationBottleneck-SteadyStateRates-TOHTOHA-",
                            "InformationBottleneck-SteadyStateRates-TOHOHA-");
      # Get the appropriate quartile from the data
      quartiles = quantile(ssr1)
      indices = switch(quartileFlag,"1"=which(ssr1 < quartiles[2]),"2"=which((ssr1 < quartiles[3])&(ssr1 >= quartiles[2])),
                       "3"=which((ssr1 < quartiles[4])&(ssr1 >= quartiles[3])),"4"=which(ssr1 >= quartiles[4]))
      ssr1 = ssr1[indices]
      ssr2 = ssr2[indices]
      
      deltaSSR = ssr2 - ssr1
      
      # Determine appropriate ROPEs to use for each parameter
      muROPE = switch(trialShortNames[1],"OH"=c(-0.3,0.3),"TOH"=c(-0.475,0.475))
      sigmaCOMP = sqrt(sd(ssr1)^2 + sd(ssr2)^2)
      sigmaROPE = c(0.9*sigmaCOMP,1.1*sigmaCOMP)
      esROPE = c(-1,1)
      
      # How precise can we be about the parameter theta, the reduction in steady-state
      # rate that occurs when a subject moves from one task to another?
      # From Jags-Ymet-Xnom1grp-Mnormal-Example.R 
      #------------------------------------------------------------------------------- 
      #------------------------------------------------------------------------------- 
      # Generate the MCMC chain:
      mcmcCoda = genMCMC( data=deltaSSR , numSavedSteps=20000 , saveName=fileNameRoot )
      #------------------------------------------------------------------------------- 
      # Display diagnostics of chain, for specified parameters:
      parameterNames = varnames(mcmcCoda) # get all parameter names
      for ( parName in parameterNames ) {
        diagMCMC( codaObject=mcmcCoda , parName=parName , 
                  saveName=fileNameRoot , saveType=graphFileType )
      }
      #------------------------------------------------------------------------------- 
      # Get summary statistics of chain:
      summaryInfo = smryMCMC( mcmcCoda , 
                              compValMu=0 , ropeMu=muROPE ,
                              compValSigma=1 , ropeSigma=sigmaROPE ,
                              compValEff=0.0 , ropeEff=esROPE ,
                              saveName=fileNameRoot )
      show(summaryInfo)
      write.csv(summaryInfo,file=paste0(fileNameRoot,"SummaryInfo-",quartileFlag,".csv"))
      # Display posterior information:
      plotMCMC( mcmcCoda , data=deltaSSR , 
                compValMu=0 , ropeMu=muROPE ,
                compValSigma=sigmaCOMP , ropeSigma=sigmaROPE ,
                compValEff=0.0 , ropeEff=esROPE ,
                pairsPlot=TRUE , showCurve=FALSE ,
                saveName=fileNameRoot , saveType=graphFileType )
      #------------------------------------------------------------------------------- 
      
      # Create a graph overlaying the three probability density functions of interest.
      healthyControlMean = switch(trialShortNames[1],"OH"=0.0665,"TOH"=0.0665)
      healthyControlSD = switch(trialShortNames[1],"OH"=0.2716,"TOH"=0.10)
      
      set.seed(1234)
      numSamples = 100000
      df <- data.frame(
        pdf=factor(rep(c("1-Region of Practical Equivalence","2-Healthy Control Model","3-Posterior Distribution"), each=numSamples)),
        samples=c(runif(numSamples, min=-muROPE[2], max=muROPE[2]),
                  rnorm(numSamples, mean=healthyControlMean, sd=healthyControlSD),
                  rnorm(numSamples, mean=summaryInfo[1,3], sd=summaryInfo[2,3]))
      )
      head(df)
      
      means <- ddply(df, "pdf", summarise, grp.mean=mean(samples))
      head(means)
      
      cols <- c("1-Region of Practical Equivalence" = "#58d68d", "2-Healthy Control Model" = "#495ab0", "3-Posterior Distribution" ="#e84133")
      
      ### GRAPH ONE
      # Use semi-transparent fill
      openGraph(height=10,width=15)
      p1<-ggplot(df[df$pdf %in% "1-Region of Practical Equivalence",], aes(x=samples, fill=pdf)) +
       coord_cartesian(xlim = c(-1.5, 1.5), ylim = c(0.0, 1.8)) +
       geom_area(aes(y = ..density..),stat ="bin", position="identity", alpha=0.4) +
       theme(legend.title = element_blank()) +
       labs(title = "Comparing Density Functions",
            subtitle = paste(trialLongNames[1]," vs. ",trialLongNames[2]),
            x = "Difference in Steady State Rates (Hz)", y = "Probability Density",
            colour = "Probability Density Function", fill = "Probability Density Function") + 
       scale_colour_manual(values = cols,aesthetics = c("colour", "fill"))
      p1
      
      # Add mean lines
      p1+geom_vline(data=means[means$pdf %in% "1-Region of Practical Equivalence",], aes(xintercept=grp.mean, color=pdf),
                   linetype="dashed")
      
      saveGraph(paste0(fileNameRoot,"PDFs-1"),graphFileType)
      
      ### GRAPH TWO
      # Use semi-transparent fill
      openGraph(height=10,width=15)
      p2<-ggplot(df[df$pdf %in% c("1-Region of Practical Equivalence","2-Healthy Control Model"),], aes(x=samples, fill=pdf)) +
       coord_cartesian(xlim = c(-1.5, 1.5), ylim = c(0.0, 1.8)) +
       geom_area(aes(y = ..density..),stat ="bin", position="identity", alpha=0.4) +
       theme(legend.title = element_blank()) +
       labs(title = "Comparing Density Functions",
            subtitle = paste(trialLongNames[1]," vs. ",trialLongNames[2]),
            x = "Difference in Steady State Rates (Hz)", y = "Probability Density",
            colour = "Probability Density Function", fill = "Probability Density Function") + 
       scale_colour_manual(values = cols,aesthetics = c("colour", "fill"))
      p2
      
      # Add mean lines
      p2+geom_vline(data=means[means$pdf %in% c("1-Region of Practical Equivalence","2-Healthy Control Model"),], aes(xintercept=grp.mean, color=pdf),
                  linetype="dashed")
      
      saveGraph(paste0(fileNameRoot,"PDFs-2"),graphFileType)
      
      ### GRAPH THREE
      # Use semi-transparent fill
      openGraph(height=10,width=15)
      
      p3<-ggplot(df, aes(x=samples, fill=pdf)) +
        coord_cartesian(xlim = c(-1.5, 1.5), ylim = c(0.0, 1.8)) +
        geom_area(aes(y = ..density..),stat ="bin", position="identity", alpha=0.6) +
        theme(legend.position = "none", axis.text = element_text(size = 72), axis.line = element_line(size = 6)) +
        labs(x="",y="") + 
        scale_colour_manual(values = cols,aesthetics = c("colour", "fill"))
      
      show(p3)
      saveGraph(paste0(fileNameRoot,"PDFs-3-",quartileFlag),graphFileType)
    }
  }
}