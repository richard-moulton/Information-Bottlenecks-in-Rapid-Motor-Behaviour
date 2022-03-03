%% Processes a single OH/OHA/TOH/TOHA trial
%  This script converts a single trial into an event/state-based format,
%  records known decisions, sets aside relevant metadata, and saves all of 
%  these for later analysis.
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
%  Requires: Communications Toolbox

function [decisionList,motorMetaData,steadyStateRate,overwhelmedHitRateMedian,overwhelmedHitRateIQR,steadyStateFrame,atLimitFrame,...
    decisionFrames,avgGap,targetCreates,distractorCreates,hits,misses,drops,crashes,subjectAccuracy,earlyBinAccuracy,...
    overwhelmedBinAccuracy,earlyBinAvoidance,overwhelmedBinAvoidance] = processSingleOHFile(data,graphFlag,~,verbose)

% We either use previous analysis to fill in needed metadata or we can make
% educated guesses based on the task type.
if isfield(data,'analysis') && isfield(data.analysis,'REPORT_FEATURES')
    motorMetaData = data.analysis.REPORT_FEATURES;
    numTargets = motorMetaData.TOTAL_TARGETS;
    
    if contains(data.file_label,["avoid","Avoid"])
        numObjects = motorMetaData.TOTAL_TARGETS + motorMetaData.TOTAL_DISTRACTORS;
    else
        numObjects = motorMetaData.TOTAL_TARGETS;
    end
else
    if verbose
        fprintf("No data.analysis or no data.analysis.REPORT_FEATURES.\n");
    end
    
    motorMetaData = [];
    
    % Hard-coded values to make up for the fact that no analysis has been done on this file.
    numObjects = 300;
    if contains(data.file_label,["avoid","Avoid"])
        numTargets = 200;
    else
        numTargets = 300;
    end
end

% Create event/state-based representations of the trial and records the
% decisions made by the subject
[eventsList,objectsList] = getEventsList(data,numObjects);
[targetCreates,distractorCreates,hits,misses,drops,crashes] = getBasicStatistics(data,eventsList);
[steadyStateRate,overwhelmedHitRateMedian,overwhelmedHitRateIQR,steadyStateFrame,atLimitFrame,subjectAccuracy] = getSteadyStateInfo(targetCreates,distractorCreates,hits,misses,drops,crashes,graphFlag,verbose);
[earlyBinAccuracy,overwhelmedBinAccuracy,earlyBinAvoidance,overwhelmedBinAvoidance] = getBinAccuracy(eventsList,objectsList,steadyStateFrame,verbose);
statesList = getStatesList(eventsList,objectsList,data);
[decisionList,decisionFrames,gaps] = getDecisionList(eventsList,statesList,atLimitFrame);

avgGap = mean(gaps);                        % Calculate the average time gap between state and decision

if verbose
    fprintf("This trial saw %i/%i targets hit (+%i/-%i), a baseline accuracy rate of %0.1f (%0.1f/%i targets) and a steady state rate of %f achieved. The average gap between state and action was %f (min %f, max %f).\n",...
        sum(hits),numTargets,sum(hits)+sum(drops),sum(misses)+sum(crashes),subjectAccuracy,...
        subjectAccuracy*numTargets,numTargets,steadyStateRate,avgGap,min(gaps),max(gaps));
end
end
