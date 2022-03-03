%% Characterize all trials by their respective task type.
%  The user can decide whether to analyse all trials, all pairs of trials,
%  or the first pair of trials.
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

% Script-level settings
verbose = 1;            % Flag, if 1 print to screen through out, if 0 don't.
graphFlag = 1;          % Flag, if 1 produce figures through out, if 0 don't.
singleTrialFlag = 2;    % Flag, set to 2 because we are analysing all trials, not just paired trials
fileSaveType = "png";   % Format in which to save figures

% Gain access to utility functions
addpath('../Utilities/');

% Decide which trials to analyse
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

% Initialize the names of variables and tasks
varNames = {'Mean hand speed left','Mean hand speed right','Mean hand speed bias','Movement area left','Movement area right','Movement area bias',...
    'Hand bias of hits','Hand selection overlap','Hand transition','Miss bias','Distractor proportion','Median error','Object processing rate',...
    'SSR (Linear)','SSR (Quadratic)'};
taskNames = {'OH','OHA','TOH','TOHA'};
r2Mat = [];

%% Object-Hit
taskName = 'Object-Hit';
taskShortName = 'OH';

temp1 = load('OHCombinedValues.mat','hitsScores','dropsScores','steadyStateRates','metaDataArray');
temp2 = load('decisionListsCompOH.mat','hitsScores','dropsScores','steadyStateRates','metaDataArray');
data = mergeStructs(temp1,temp2);

if graphFlag
    f = hitsScoresHistogram(data.hitsScores,taskName); %#ok<*UNRCH>
    fileNameString = sprintf('rangeOfHits(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    
    lmdl = fitlm(data.steadyStateRates,data.hitsScores);
    f = blandAltmanPlot(data.hitsScores,predict(lmdl,data.steadyStateRates),...
        'Average of Targets Hit','Difference between Targets Hit','');
    fileNameString = sprintf('blandAltmanLinear(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    fprintf('Linear model for %s: NumHits = %f SSR + %f\n',taskName,lmdl.Coefficients.Estimate(2),lmdl.Coefficients.Estimate(1));
    
    qmdl = fitlm(data.steadyStateRates,data.hitsScores,'quadratic');
    f = blandAltmanPlot(data.hitsScores,predict(qmdl,data.steadyStateRates),...
        'Average of Targets Hit','Difference between Targets Hit','');
    fileNameString = sprintf('blandAltmanQuadratic(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    fprintf('Quadratic model for %s: NumHits = %f SSR^2 + %f SSR + %f\n',taskName,qmdl.Coefficients.Estimate(3),qmdl.Coefficients.Estimate(2),qmdl.Coefficients.Estimate(1));
end

r2 = characterizeSingleTaskType(data,fileSaveType,taskName,taskShortName,singleTrialFlag,verbose,graphFlag);
r2Mat = [r2Mat;r2];

close all;
pause(1);

%% Object-Hit-and-Avoid
taskName = 'Object-Hit-and-Avoid';
taskShortName = 'OHA';

data = load('OHACombinedValues.mat');

if graphFlag
    f = hitsScoresHistogram(data.hitsScores,taskName);
    fileNameString = sprintf('rangeOfHits(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    
    lmdl = fitlm(data.steadyStateRates,data.hitsScores);
    f = blandAltmanPlot(data.hitsScores,predict(lmdl,data.steadyStateRates),...
        'Average of Targets Hit','Difference between Targets Hit','');
    fileNameString = sprintf('blandAltmanLinear(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    fprintf('Linear model for %s: NumHits = %f SSR + %f\n',taskName,lmdl.Coefficients.Estimate(2),lmdl.Coefficients.Estimate(1));
    
    qmdl = fitlm(data.steadyStateRates,data.hitsScores,'quadratic');
    f = blandAltmanPlot(data.hitsScores,predict(qmdl,data.steadyStateRates),...
        'Average of Targets Hit','Difference between Targets Hit','');
    fileNameString = sprintf('blandAltmanQuadratic(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    fprintf('Quadratic model for %s: NumHits = %f SSR^2 + %f SSR + %f\n',taskName,qmdl.Coefficients.Estimate(3),qmdl.Coefficients.Estimate(2),qmdl.Coefficients.Estimate(1));
end

r2 = characterizeSingleTaskType(data,fileSaveType,taskName,taskShortName,singleTrialFlag,verbose,graphFlag);
r2Mat = [r2Mat;r2];

close all;
pause(1);

%% Turbo Object-Hit
taskName = 'Turbo Object-Hit';
taskShortName = 'TOH';

temp1 = load('turboOHCombinedValues.mat');
temp2 = load('decisionListsCompTOH.mat');
data = mergeStructs(temp1,temp2);

if graphFlag
    f = hitsScoresHistogram(data.hitsScores,taskName);
    fileNameString = sprintf('rangeOfHits(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    
    lmdl = fitlm(data.steadyStateRates,data.hitsScores);
    f = blandAltmanPlot(data.hitsScores,predict(lmdl,data.steadyStateRates),...
        'Average of Targets Hit','Difference between Targets Hit','');
    fileNameString = sprintf('blandAltmanLinear(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    fprintf('Linear model for %s: NumHits = %f SSR + %f\n',taskName,lmdl.Coefficients.Estimate(2),lmdl.Coefficients.Estimate(1));
    
    qmdl = fitlm(data.steadyStateRates,data.hitsScores,'quadratic');
    f = blandAltmanPlot(data.hitsScores,predict(qmdl,data.steadyStateRates),...
        'Average of Targets Hit','Difference between Targets Hit','');
    fileNameString = sprintf('blandAltmanQuadratic(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    fprintf('Quadratic model for %s: NumHits = %f SSR^2 + %f SSR + %f\n',taskName,qmdl.Coefficients.Estimate(3),qmdl.Coefficients.Estimate(2),qmdl.Coefficients.Estimate(1));
end

r2 = characterizeSingleTaskType(data,fileSaveType,taskName,taskShortName,singleTrialFlag,verbose,graphFlag);
r2Mat = [r2Mat;r2];

close all;
pause(1);

%% Turbo Object-Hit-and-Avoid
taskName = 'Turbo Object-Hit-and-Avoid';
taskShortName = 'TOHA';

data = load('turboOHACombinedValues.mat');

if graphFlag
    f = hitsScoresHistogram(data.hitsScores,taskName);
    fileNameString = sprintf('rangeOfHits(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    
    lmdl = fitlm(data.steadyStateRates,data.hitsScores);
    f = blandAltmanPlot(data.hitsScores,predict(lmdl,data.steadyStateRates),...
        'Average of Targets Hit','Difference between Targets Hit','');
    fileNameString = sprintf('blandAltmanLinear(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    fprintf('Linear model for %s: NumHits = %f SSR + %f\n',taskName,lmdl.Coefficients.Estimate(2),lmdl.Coefficients.Estimate(1));
    
    qmdl = fitlm(data.steadyStateRates,data.hitsScores,'quadratic');
    f = blandAltmanPlot(data.hitsScores,predict(qmdl,data.steadyStateRates),...
        'Average of Targets Hit','Difference between Targets Hit','');
    fileNameString = sprintf('blandAltmanQuadratic(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    fprintf('Quadratic model for %s: NumHits = %f SSR^2 + %f SSR + %f\n',taskName,qmdl.Coefficients.Estimate(3),qmdl.Coefficients.Estimate(2),qmdl.Coefficients.Estimate(1));
end

r2 = characterizeSingleTaskType(data,fileSaveType,taskName,taskShortName,singleTrialFlag,verbose,graphFlag);
r2Mat = [r2Mat;r2];

close all;
pause(1);

if graphFlag
    f = figure;
    h = heatmap(taskNames,varNames,r2Mat','Colormap',jet,'ColorLimits',[0 1]);
    h.Title = 'Adjusted R^{2} values for subject variables and number of targets hit';
    h.XLabel = 'Task';
    h.YLabel = 'Variables';
    ax = gca;
    ax.FontSize = 24;
    f.WindowState = 'maximized';
    fileNameString = sprintf('r2heatMap%s',singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    
    f = figure;
    h = heatmap(taskNames,varNames(1:10),r2Mat(:,1:10)','Colormap',jet,'ColorLimits',[0 1]);
    h.Title = 'Adjusted R^{2} values for subject variables and number of targets hit';
    h.XLabel = 'Task';
    h.YLabel = 'Variables';
    ax = gca;
    ax.FontSize = 24;
    f.WindowState = 'maximized';
    fileNameString = sprintf('r2heatMapMotorControl%s',singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
end

function f = hitsScoresHistogram(hitsScores,taskName)
switch taskName
    case 'Object-Hit'
        maxHits = 300;
    case 'Object-Hit-and-Avoid'
        maxHits = 200;
    case 'Turbo Object-Hit'
        maxHits = 300;
    case 'Turbo Object-Hit-and-Avoid'
        maxHits = 200;
    otherwise
        fprintf('Unexpected task name! Got ''%s''.\n',taskName);
        keyboard
end

hitsRange = 0:5:maxHits;

meanScore = mean(hitsScores);
stdScore = std(hitsScores);
medianScore = median(hitsScores);

fprintf('%s\nThe mean subject hit %f targets. The standard deviation across the data set was %f.\n',taskName,meanScore,stdScore);
fprintf('The median subject hit %f targets.\n',medianScore);

f = figure;
hold on;
dataString = sprintf('Data (n = %i)',length(hitsScores));
histogram(hitsScores,hitsRange,'FaceColor',[0 0.4470 0.7410],'FaceAlpha',0.3,'DisplayName',dataString)
xline(maxHits,'-k','LineWidth',4,'DisplayName','Perfect Performance');
xline(meanScore,'-','Color',[1 0.2 0.2],'LineWidth',4,'DisplayName','Mean');
xline(meanScore + stdScore,'--','Color',[1 0.2 0.2],'LineWidth',2,'DisplayName','+/- Std Dev');
xline(meanScore - stdScore,'--','Color',[1 0.2 0.2],'LineWidth',2,'HandleVisibility','off');
xline(medianScore,'-','Color',[0 0.6 0.3],'LineWidth',4,'DisplayName','Median');
lgd = legend('location','northwest');
ylabel("Number of Subjects (n)");
xlabel("Number of Targets Hit (n)");
axis([0 maxHits 0 Inf]);
ax = gca;
ax.FontSize = 32;
ax.FontWeight = 'bold';
ax.LineWidth = 2.5;
lgd.FontSize = 24;
lgd.FontWeight = 'normal';
f.WindowState = 'maximized';
hold off;
end
