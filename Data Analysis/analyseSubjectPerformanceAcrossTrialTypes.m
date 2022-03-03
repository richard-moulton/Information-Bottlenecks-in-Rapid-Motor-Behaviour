%% Compares performance in two tasks on a "per-subject" basis
%  Given two tasks, this script compares them by analysing the performance
%  of subjects who performed both tasks. The user has the option of
%  comparing this performance between all aligned pairs of trials for
%  subjects, or only each subject's first pair of trials.
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
%  8 December, 2021
%
%  Requires: Statistics and Machine Learning Toolbox
%#ok<*UNRCH>

% Script-level parameters
compFlag = 4;           % Flag indicating which comparison to make
verbose = 1;            % Flag, if 1 print to screen through out, if 0 don't.
graphFlag = 1;          % Flag, if 1 produce figures through out, if 0 don't.
singleTrialFlag = 1;    % Flag, if 1 take only the first trial pair for each subject, if 0 take all aligned trial pairs
fileSaveType = "png";   % Format in which to save figures

% Gain access to utility functions
addpath('../Utilities/');

%% Load data and set variables
switch compFlag
    case 1
        t1Results = load('OHCombinedValues.mat');
        t2Results = load('OHACombinedValues.mat');
        dayMatchFlag = 1;
    case 2
        t1Results = load('decisionListsCompOH.mat');
        t2Results = load('decisionListsCompTOH.mat');
        dayMatchFlag = 0;
    case 3
        t1Results = load('turboOHCombinedValues.mat');
        t2Results = load('turboOHACombinedValues.mat');
        dayMatchFlag = 1;
    case 4
        t1Results = load('decisionListsCompTOH.mat');
        t2Results = load('OHACombinedValues.mat');
        dayMatchFlag = 0;
    otherwise
        fprintf('Bad value for compFlag!');
        keyboard
end

% Set task names
if t1Results.OHflag
    t1Name = 'Object-Hit';
    t1ShortName = 'OH';
else
    t1Name = 'Object-Hit-and-Avoid';
    t1ShortName = 'OHA';
end
if t2Results.OHflag
    t2Name = 'Object-Hit';
    t2ShortName = 'OH';
else
    t2Name = 'Object-Hit-and-Avoid';
    t2ShortName = 'OHA';
end
if t1Results.turboFlag
    t1Name = ['Turbo ' t1Name];
    t1ShortName = ['T' t1ShortName];
end
if t2Results.turboFlag
    t2Name = ['Turbo ' t2Name];
    t2ShortName = ['T' t2ShortName];
end

% Set string modifier for single trials
if singleTrialFlag
    singleTrialString = "-SingleTrial";
else
    singleTrialString = "";
end

%% Determine which subjects in trial 1 have also performed trial 2.
%  From there, determine based on trial time-stamps how to align the
%  trials.
t1UniqueSubjects = unique(t1Results.subjectNames);

estLength = length(t1Results.steadyStateRates);
trial1SSRs = zeros(estLength,1);
trial2SSRs = zeros(estLength,1);
trial1Max = zeros(estLength,1);
trial2Max = zeros(estLength,1);
idx = 1;
numSubjectsMatched = 0;
subjectsMatched = strings(estLength,1);

trial1SDs = NaN(length(t1UniqueSubjects),1);

for i=1:length(t1UniqueSubjects)
    subjectName = t1UniqueSubjects(i);
    
    if verbose
        fprintf('Looking for subject %s...',subjectName);
    end
    
    t1TrialIndices = find(ismember(t1Results.subjectNames,subjectName));
    t2TrialIndices = find(ismember(t2Results.subjectNames,subjectName));
    
    if length(t1TrialIndices) > 1
        trial1SDs(i) = std(t1Results.steadyStateRates(t1TrialIndices));
    end
    
    [t1TrialIndices,t2TrialIndices] = alignTrials(t1TrialIndices,t2TrialIndices,t1Results.trialNames,t2Results.trialNames,singleTrialFlag,dayMatchFlag,verbose);
    
    if ~isempty(t1TrialIndices) && ~isempty(t2TrialIndices)
        numSubjectsMatched = numSubjectsMatched + 1;
    end
    
    newIdx = idx + length(t1TrialIndices) - 1;
    trial1SSRs(idx:newIdx) = t1Results.steadyStateRates(t1TrialIndices);
    trial2SSRs(idx:newIdx) = t2Results.steadyStateRates(t2TrialIndices);
    trial1Max(idx:newIdx) = t1Results.maxAvgTargetCreationRate(t1TrialIndices);
    trial2Max(idx:newIdx) = t2Results.maxAvgTargetCreationRate(t2TrialIndices);
    subjectsMatched(idx:newIdx) = subjectName;
    idx = newIdx + 1;
    
    if verbose
        fprintf("There are %i matches in Trial 1 and %i matches in Trial 2.\n",length(t1TrialIndices),length(t2TrialIndices));
        for j=1:length(t1TrialIndices)
            fprintf("%s - %s\n",t1Results.trialNames(t1TrialIndices(j)), t2Results.trialNames(t2TrialIndices(j)));
        end
    end
end

% Trim arrays as necessary
trial1SSRs = trial1SSRs(1:idx-1);
trial2SSRs = trial2SSRs(1:idx-1);
trial1Max = trial1Max(1:idx-1);
trial2Max = trial2Max(1:idx-1);
subjectsMatched = subjectsMatched(1:idx-1);

% Get the peak target creation rate across both trial types
maxTargetCreationRate = max([trial1Max;trial2Max]);

% Get the minimum SSRs across both trial types
minSSR = min([trial1SSRs;trial2SSRs]);
maxSSR = max([trial1SSRs;trial2SSRs]);

% Get the mean standard deviation of SSRs for subjects with more than one
% t1 trial
if sum(~isnan(trial1SDs)) > 10
    t1MeanSD = mean(trial1SDs,'omitnan');
    t1ROPE = 2*t1MeanSD;
elseif strcmp(t1ShortName,'OH')
    t1MeanSD = 0.15;
    t1ROPE = 0.30;
elseif strcmp(t1ShortName,'TOH')
    t1MeanSD = 0.24;
    t1ROPE = 0.48;
end

% If there are subjects who have performed these tasks more than once,
% investigate whether a subject's variability in performance is related to
% their skill level.
x = trial1SSRs(~isnan(trial1SDs));
y = trial1SDs(~isnan(trial1SDs));
if ~isempty(x) && ~isempty(y)
    [f,mdl] = scatterWithRegression(x,y,[0 4],[0 1],1);
    set(groot,'currentfigure',f);
    hold on;
    hold off;
    fprintf('The adjusted R2 value for SSR standard deviations vs SSR for %s is %0.2f. The slope is %0.2f.\n',t1Name,mdl.Rsquared.Adjusted,mdl.Coefficients.Estimate(2));
end

% Produce useful strings using the trial task types
trialSubjectString = sprintf('(%i trials, %i subjects)',length(trial1SSRs),numSubjectsMatched);
maxString1 = sprintf('Peak Target Creation Rates (%s)',t1ShortName);
maxString2 = sprintf('Peak Target Creation Rates (%s)',t2ShortName);
ropeString = sprintf('The mean standard deviation for subjects with more than one %s trial is %1.5f. This gives a ROPE of %1.5f.',t1Name,t1MeanSD,t1ROPE);
fprintf('%s\n',ropeString);

%% Boxplot and scatterplot for raw steady-state rates
if graphFlag
    group = [    ones(size(trial1SSRs,1),1);
        2 * ones(size(trial2SSRs,1),1)];
    
    f = figure;
    hold on;
    titleString = sprintf('Comparing within-subject steady-state rates\nin %s and %s %s',...
        t1Name,t2Name,trialSubjectString);
    title(titleString);
    a = yline(max(trial1Max),'--r','LineWidth',2);
    b = yline(min(trial1Max),'--r','LineWidth',2);
    c = yline(max(trial2Max),':r','LineWidth',2);
    d = yline(min(trial2Max),':r','LineWidth',2);
    boxplot([trial1SSRs; trial2SSRs],group,'labels',{t1Name,t2Name});
    e = parallelcoords([trial1SSRs, trial2SSRs], 'Color', 0.7*[1 1 1], 'LineStyle', '-',...
        'Marker', '.', 'MarkerSize', 10);
    ylim([0 ceil(maxTargetCreationRate)]);
    legend([a(1) c(1) e(1)],maxString1,maxString2,'Paired Observations');
    xlabel('Task');
    ylabel('Steady state rates (Hz)');
    ax = gca;
    ax.FontSize = 24;
    f.WindowState = 'maximized';
    hold off;
    fileNameString = sprintf('intrasubjectSteadyStateRate(%s-%s%s)',t1ShortName,t2ShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
end

f = figure;
hold on;
scatter(trial1SSRs,trial2SSRs,'filled','b','DisplayName',sprintf('Participants (%i)',length(trial1SSRs)));
plot([minSSR maxSSR],[minSSR maxSSR],'-k','LineWidth',2,'DisplayName','Unity Line');
plot([minSSR+t1ROPE maxSSR+t1ROPE],[minSSR maxSSR],'-b','DisplayName','ROPE');
plot(minSSR:0.1:maxSSR,predict(fitlm(trial1SSRs,trial2SSRs),(minSSR:0.1:maxSSR)'),'-r','LineWidth',2,'DisplayName','Regression');
lgd = legend('location','best');
xlabel(sprintf('%s Steady State Rate (Hz)',t1ShortName));
ylabel(sprintf('%s Steady State Rate (Hz)',t2ShortName));
ax = gca;
ax.FontSize = 32;
ax.FontWeight = 'bold';
lgd.FontSize = 24;
lgd.FontWeight = 'normal';
f.WindowState = 'maximized';
hold off;
fileNameString = sprintf('intrasubjectSteadyStateRateScatter(%s-%s%s)',t1ShortName,t2ShortName,singleTrialString);
savefigas(f,fileNameString,fileSaveType);

if verbose
    fprintf("Mean difference between steady-state rates in %s and %s: %0.3f\n",t1ShortName,t2ShortName,mean(trial1SSRs - trial2SSRs));
end

%% Boxplot for adjusted steady-state rates
adjustedTrial1SSRs = trial1SSRs.*maxTargetCreationRate./trial1Max;
adjustedTrial2SSRs = trial2SSRs.*maxTargetCreationRate./trial2Max;

if graphFlag
    f = figure;
    hold on;
    titleString = sprintf('Comparing within-subject steady-state rates adjusted for\ntarget creation rates in %s and %s %s',...
        t1Name,t2Name,trialSubjectString);
    title(titleString);
    yline(maxTargetCreationRate,'r','LineWidth',2);
    boxplot([adjustedTrial1SSRs; adjustedTrial2SSRs],...
        group,'labels',{t1Name,t2Name});
    parallelcoords([adjustedTrial1SSRs, adjustedTrial2SSRs], 'Color', 0.7*[1 1 1], 'LineStyle', '-',...
        'Marker', '.', 'MarkerSize', 10);
    ylim([0 ceil(maxTargetCreationRate)]);
    legend('Peak Target Creation Rate','Paired Observations');
    xlabel('Task');
    ylabel('Adjusted Steady state rates (Hz)');
    ax = gca;
    ax.FontSize = 24;
    f.WindowState = 'maximized';
    hold off;
    fileNameString = sprintf('intrasubjectAdjustedSteadyStateRate(%s-%s%s)',t1ShortName,t2ShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
end

if verbose
    fprintf("Mean difference between adjusted steady-state rates in %s and %s: %0.3f\n",t1ShortName,t2ShortName,mean(trial1SSRs.*maxTargetCreationRate./trial1Max - trial2SSRs.*maxTargetCreationRate./trial2Max));
end

f = figure;
hold on;
scatter(adjustedTrial1SSRs, adjustedTrial2SSRs,'filled','b','DisplayName','Male');
plot([minSSR maxSSR],[minSSR maxSSR],'-k','LineWidth',2,'DisplayName','Unity Line');
plot([minSSR+t1ROPE maxSSR+t1ROPE],[minSSR maxSSR],'-b','DisplayName','ROPE');
plot(minSSR:0.1:maxSSR,predict(fitlm(adjustedTrial1SSRs, adjustedTrial2SSRs),(minSSR:0.1:maxSSR)'),'-r','LineWidth',2,'DisplayName','Regression');
lgd = legend('location','best');
xlabel(sprintf('%s Adjusted Steady State Rate (Hz)',t1ShortName));
ylabel(sprintf('%s Adjusted Steady State Rate (Hz)',t2ShortName));
ax = gca;
ax.FontSize = 32;
ax.FontWeight = 'bold';
lgd.FontSize = 24;
lgd.FontWeight = 'normal';
f.WindowState = 'maximized';
hold off;
fileNameString = sprintf('intrasubjectAdjustedSteadyStateRateScatter(%s-%s%s)',t1ShortName,t2ShortName,singleTrialString);
savefigas(f,fileNameString,fileSaveType);


%% Scatterplot relating steady-state rates between tasks
if graphFlag
    titleString = sprintf('Steady state rates correlate between %s and %s',t1Name,t2Name);
    t1Label = sprintf('%s Steady State Rate (Hz)',t1ShortName);
    t2Label = sprintf('%s Steady State Rate (Hz)',t2ShortName);
    
    f = scatterAndRegressionWithCIs(trial1SSRs,trial2SSRs,length(trial1SSRs),...
        t1Label,t2Label,titleString,'');
    hold on;
    plot([minSSR maxSSR],[minSSR maxSSR],'-k','DisplayName',"Unity Line");
    hold off;
    fileNameString = sprintf('SteadyStateRateCorrelations(%s-%s%s)',t1ShortName,t2ShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
end

%% Save steady-state rates as a csv file for Bayesian Data Analysis in R
fileNameString = sprintf('SubjectSteadyStateRates(%s-%s%s).csv',t1ShortName,t2ShortName,singleTrialString);
csvwrite(fileNameString,[trial1SSRs trial2SSRs]);

close all;

function [f,mdl] = scatterWithRegression(x,y,xLimits,yLimits,unityFlag)
f = figure;
hold on;

% Scatter
scatter(x,y,150,'black','filled','MarkerFaceAlpha',0.6);
pause(1);

% Unity Line
if unityFlag
    minCoord = min([xLimits(1)  yLimits(1)]);
    maxCoord = max([xLimits(2) yLimits(2)]);
    plot([minCoord maxCoord],[minCoord maxCoord],'k','LineWidth',6)
else
    plot([0 200],[0 100],'k','LineWidth',6);
end

% Regression
mdl = fitlm(x,y,'linear');
xVals = linspace(min(x),max(x));
[mdlVals,~] = predict(mdl,xVals','Alpha',0.05);

plot(xVals,mdlVals,'Color',[0.8500 0.3250 0.0980],'LineWidth',10);

%title(sprintf('n = %i',length(x)));

xlim(xLimits);
ylim(yLimits);

ax = gca;
ax.FontSize = 72;
ax.FontWeight = 'bold';
ax.LineWidth = 15;
axis square
f.WindowState = 'maximized';
hold off;
end
