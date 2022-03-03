%% Characterize performance in a single task
%  The question answered by this function is, which recorded variable best
%  explains task performance (i.e., the number of targets hit by the
%  participant)?
%  Copyright (C) 2022  Richard Hugh Moulton
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%  6 December, 2021
%
%  Requires: Statistics and Machine Learning Toolbox

function r2Values = characterizeSingleTaskType(data,fileSaveType,taskName,taskShortName,singleTrialFlag,verbose,graphFlag)

% singleTrialFlag is overloaded to denote the kind of trials we have.
switch singleTrialFlag
    case 2  % All trials, not just pairs
        singleTrialString = '-AllTrials';
    case 1  % Trials from all first pairs
        singleTrialString = '-SingleTrial';
    case 0  % Trials from all pairs
        singleTrialString = '-AllPairs';
    otherwise
        fprintf('Bad value for singleTrialFlag! Got %i and was expecting 0, 1, or 2.\n',singleTrialFlag);
end

if verbose
    fprintf("Characterizing %s\n",taskName);
end

% Initialize arrays and gain access to utility functions
addpath('../Utilities/');
r2Values = NaN(1,15);

% Calculate some global variables relating to all trials
muSSR = mean(data.steadyStateRates);
sdSSR = std(data.steadyStateRates);
numTrials = length(data.steadyStateRates);

N = size(data.steadyStateRates,1);                  % Number of ‘Experiments’ In Data Set
ssrSEM = sdSSR/sqrt(N);                             % Compute ‘Standard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
ssrCI95 = bsxfun(@times, ssrSEM, CI95(:));  

normalizedSSR = (data.steadyStateRates - muSSR)./sdSSR;
[h,p] = kstest(normalizedSSR);

if verbose
    fprintf('Steady state rates over %i %s trials have a mean of %0.6f and a standard deviation of %0.6f.\n',numTrials,taskName,muSSR,sdSSR);
    
    fprintf('The 95%% confidence interval for the mean steady state rate is [%0.6f, %0.6f].\n',ssrCI95);
    
    if h
        fprintf('The Kolmogorov-Smirnov test rejected the null hypothesis that this data is normally distributed (p-value %0.6f).\n',p);
    else
        fprintf('The Kolmogorov-Smirnov test failed to reject the null hypothesis that this data is normally distributed (p-value %0.6f).\n',p);
    end
end

distractors = data.metaDataArray.HAS_DISTRACTORS;
if graphFlag && strcmp(distractors,'true')
    % Do hits and drops correlate?
    titleString = sprintf('Subjects who hit targets also avoid distractors in %s',taskName);
    xLabel = sprintf('Number of Targets Hit (n)');
    yLabel = sprintf('Number of Distractors Avoided (n)');
    
    f = scatterAndRegressionWithCIs(data.hitsScores,data.dropsScores,length(data.hitsScores),...
        xLabel,yLabel,titleString,'');
    hold on;
    plot([0 200],[0 100],'-k','DisplayName',"Equal Performance");
    hold off;
    fileNameString = sprintf('HitsDropsCorrelations(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
end

% Need to find trials where we have meta data array!
m = data.metaDataArray;
fullIdx = find(arrayfun(@(m) ~isempty(m.DISTRACTOR_HITS_LEFT),m));
hits = [data.hitsScores];
hits = hits(fullIdx);
ssrs = [data.steadyStateRates];
ssrs = ssrs(fullIdx);

%% Do motor control characteristics explain performance (total hits?)
[f,mdl] =    scatterAndRegressionWithCIs([data.metaDataArray.HAND_BIAS_OF_HITS],hits,numTrials,...
    'Hand Bias of Hits','Targets Hits (n)','Hand bias of hits does not predict overall task performance in ',taskName);
r2Values(7) = mdl.Rsquared.Adjusted;
if verbose
    fprintf("Linear regression between hits and hand bias of hits: R^{2} = %0.2f\n",mdl.Rsquared.Adjusted);
end
if graphFlag
    fileNameString = sprintf('HandBiasofHits(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
end

[f,mdl] = scatterAndRegressionWithCIs([data.metaDataArray.HAND_SELECTION_OVERLAP],hits,numTrials,...
    'Hand Selection Overlap (% of hits)','Targets Hits (n)','Hand selection overlap does not predict overall task performance in ',taskName);
r2Values(8) = mdl.Rsquared.Adjusted;
if verbose
    fprintf("Linear regression between hits and hand selection overlap: R^{2} = %0.2f\n",mdl.Rsquared.Adjusted);
end
if graphFlag
    fileNameString = sprintf('HandSelectionOverlap(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
end

[f,mdl] =    scatterAndRegressionWithCIs([data.metaDataArray.HAND_TRANSITION],hits,numTrials,...
    'Hand Transition (cm)','Targets Hits (n)','Hand transition point does not predict overall task performance in ',taskName);
r2Values(9) = mdl.Rsquared.Adjusted;
if verbose
    fprintf("Linear regression between hits and hand transition: R^{2} = %0.2f\n",mdl.Rsquared.Adjusted);
end
if graphFlag
    fileNameString = sprintf('HandTransition(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
end

[f,mdl] =    scatterAndRegressionWithCIs([data.metaDataArray.MEAN_HAND_SPEED_BIAS],hits,numTrials,...
    'Mean Hand Speed Bias','Targets Hits (n)','Mean hand speed bias does not predict overall task performance in ',taskName);
r2Values(3) = mdl.Rsquared.Adjusted;
if verbose
    fprintf("Linear regression between hits and mean hand speed bias: R^{2} = %0.2f\n",mdl.Rsquared.Adjusted);
end
if graphFlag
    fileNameString = sprintf('MeanHandSpeedBias(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
end

[f,mdl] =    scatterAndRegressionWithCIs([data.metaDataArray.MEAN_HAND_SPEED_LEFT],hits,numTrials,...
    'Mean Hand Speed Left (cm/s)','Targets Hits (n)','Mean hand speed left does not predict overall task performance in ',taskName);
r2Values(1) = mdl.Rsquared.Adjusted;
if verbose
    fprintf("Linear regression between hits and mean hand speed left: R^{2} = %0.2f\n",mdl.Rsquared.Adjusted);
end
if graphFlag
    fileNameString = sprintf('MeanHandSpeedLeft(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
end

[f,mdl] =    scatterAndRegressionWithCIs([data.metaDataArray.MEAN_HAND_SPEED_RIGHT],hits,numTrials,...
    'Mean Hand Speed Right (cm/s)','Targets Hits (n)','Mean hand speed right does not predict overall task performance in ',taskName);
r2Values(2) = mdl.Rsquared.Adjusted;
if verbose
    fprintf("Linear regression between hits and mean hand speed right: R^{2} = %0.2f\n",mdl.Rsquared.Adjusted);
end
if graphFlag
    fileNameString = sprintf('MeanHandSpeedRight(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
end

[f,mdl] =    scatterAndRegressionWithCIs([data.metaDataArray.MEDIAN_ERROR],hits,numTrials,...
    'Median Error (% of targets)','Targets Hits (n)','Median error predicts overall task performance in ',taskName);
r2Values(12) = mdl.Rsquared.Adjusted;
if verbose
    fprintf("Linear regression between hits and median error: R^{2} = %0.2f\n",mdl.Rsquared.Adjusted);
end
if graphFlag
    fileNameString = sprintf('MedianError(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
end

[f,mdl] =    scatterAndRegressionWithCIs([data.metaDataArray.MISS_BIAS],hits,numTrials,...
    'Miss Bias (cm)','Targets Hits (n)','Miss bias does not predict overall task performance in ',taskName);
r2Values(10) = mdl.Rsquared.Adjusted;
if verbose
    fprintf("Linear regression between hits and miss bias: R^{2} = %0.2f\n",mdl.Rsquared.Adjusted);
end
if graphFlag
    fileNameString = sprintf('MissBias(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
end

[f,mdl] =    scatterAndRegressionWithCIs([data.metaDataArray.MOVEMENT_AREA_BIAS],hits,numTrials,...
    'Movement Area Bias','Targets Hits (n)','Movement area bias does not predict overall task performance in ',taskName);
r2Values(6) = mdl.Rsquared.Adjusted;
if verbose
    fprintf("Linear regression between hits and movement area bias: R^{2} = %0.2f\n",mdl.Rsquared.Adjusted);
end
if graphFlag
    fileNameString = sprintf('MovementAreaBias(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
end

[f,mdl] =    scatterAndRegressionWithCIs([data.metaDataArray.MOVEMENT_AREA_LEFT],hits,numTrials,...
    'Movement Area Left (cm^{2})','Targets Hits (n)','Movement area left does not predict overall task performance in ',taskName);
r2Values(4) = mdl.Rsquared.Adjusted;
if verbose
    fprintf("Linear regression between hits and movement area left: R^{2} = %0.2f\n",mdl.Rsquared.Adjusted);
end
if graphFlag
    fileNameString = sprintf('MovementAreaLeft(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
end

[f,mdl] =    scatterAndRegressionWithCIs([data.metaDataArray.MOVEMENT_AREA_RIGHT],hits,numTrials,...
    'Movement Area Right (cm^{2})','Targets Hits (n)','Movement area right does not predict overall task performance in ',taskName);
r2Values(5) = mdl.Rsquared.Adjusted;
if verbose
    fprintf("Linear regression between hits and movement area right: R^{2} = %0.2f\n",mdl.Rsquared.Adjusted);
end
if graphFlag
    fileNameString = sprintf('MovementAreaRight(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
end

%% Do steady state rates explain performance (total hits?)
[f,mdl] =    scatterAndRegressionWithCIs(ssrs,hits,numTrials,...
    'Steady State Rate (Hz)','Targets Hits (n)','Steady state rate predicts overall task performance in ',taskName);
r2Values(14) = mdl.Rsquared.Adjusted;
if verbose
    fprintf("Linear regression between hits and steady state rate: R^{2} = %0.2f\n",mdl.Rsquared.Adjusted);
end
if graphFlag
    fileNameString = sprintf('SteadyStateRateLinear(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
end

mdl = fitlm(ssrs,hits,'quadratic');
xVals = linspace(min(ssrs),max(ssrs));
[mdlVals,mdlCIs] = predict(mdl,xVals','Alpha',0.05);
r2Values(15) = mdl.Rsquared.Adjusted;

if graphFlag
    f = figure;
    hold on;
    
    scatter(ssrs,hits,'filled');
    plot(xVals,mdlVals,'Color',[0.8500 0.3250 0.0980],'LineWidth',2);
    plot(xVals,mdlCIs(:,1),'--','Color',[0.8500 0.3250 0.0980]);
    plot(xVals,mdlCIs(:,2),'--','Color',[0.8500 0.3250 0.0980]);
    
    dataString = sprintf('Data (n = %i)',numTrials);
    lrString = sprintf('Quadratic Regression (R^{2} = %1.2f)',mdl.Rsquared.Adjusted);
    lgd = legend(dataString,lrString,'95% CI','location','northeast');
    
    xlabel('Steady State Rate (Hz)');
    ylabel('Targets Hit (n)');
    ax = gca;
    ax.FontSize = 32;
    ax.FontWeight = 'bold';
    ax.LineWidth = 2.5;
    lgd.FontSize = 24;
    lgd.FontWeight = 'normal';
    f.WindowState = 'maximized';
    hold off;
    
    fileNameString = sprintf('SteadyStateRateQuadratic(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
end

if isfield(data.metaDataArray,'DISTRACTOR_PROPORTION')
    fullIdx = find(arrayfun(@(m) ~isempty(m.DISTRACTOR_PROPORTION),m));
    hits = [data.hitsScores];
    hits = hits(fullIdx);
    [f,mdl] =    scatterAndRegressionWithCIs([data.metaDataArray.DISTRACTOR_PROPORTION],hits,numTrials,...
        'Distractor Proportion (%)','Targets Hits (n)','Distractor proportion predicts overall task performance in ',taskName);
    r2Values(11) = mdl.Rsquared.Adjusted;
    if verbose
        fprintf("Linear regression between hits and distractor proportion: R^{2} = %0.2f\n",mdl.Rsquared.Adjusted);
    end
    if graphFlag
        fileNameString = sprintf('DistractorProportion(%s%s)',taskShortName,singleTrialString);
        savefigas(f,fileNameString,fileSaveType);
    end
end

if isfield(data.metaDataArray,'OBJECT_PROCESSING_RATE')
    fullIdx = find(arrayfun(@(m) ~isempty(m.OBJECT_PROCESSING_RATE),m));
    hits = [data.hitsScores];
    hits = hits(fullIdx);
    [f,mdl] =    scatterAndRegressionWithCIs([data.metaDataArray.OBJECT_PROCESSING_RATE],hits,numTrials,...
        'Object Processing Rate (Hz)','Targets Hits (n)','Object processing rate predicts overall task performance in ',taskName);
    r2Values(13) = mdl.Rsquared.Adjusted;
    if verbose
        fprintf("Linear regression between hits and object processing rate: R^{2} = %0.2f\n",mdl.Rsquared.Adjusted);
    end
    if graphFlag
        fileNameString = sprintf('ObjectProcessingRate(%s%s)',taskShortName,singleTrialString);
        savefigas(f,fileNameString,fileSaveType);
    end
end
end
