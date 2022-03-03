%% Determine a participant's steady-state performance in OH/OHA/TOH/TOA
%  Automatically detects three phases of task performance (early/easy,
%  at-limit, and overwhelmed) and calculates some phase-related statistics.
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

function [steadyStateRate,overwhelmedHitRateMedian,overwhelmedHitRateIQR,...
    steadyStateFrame,atLimitFrame,subjectAccuracy] = getSteadyStateInfo(targetCreates,distractorCreates,hits,misses,~,~,graphFlag,verbose)

% We multiply the event-occurence arrays by 200 to transform their averages
% from 'per-frame' to 'per-second'.
targetCreates = targetCreates * 200;
distractorCreates = distractorCreates * 200;
hits = hits * 200;
objectCreates = targetCreates + distractorCreates;

% Implement Kolmogorov-Zurbenko filters: 2 iterations, window sizes 1001
% data points
avgObjectCreationRate = movmean(movmean(objectCreates,1001,'omitnan'),1001,'omitnan');
avgTargetCreationRate = movmean(movmean(targetCreates,1001,'omitnan'),1001,'omitnan');
avgHitRate = movmean(movmean(hits,1001,'omitnan'),1001,'omitnan');

% Calculate where the average object creation rate begins to fall off, this
% is the "end of the task" where the subject is no longer being
% increasingly stressed.
[~,endOfTask] = max(avgObjectCreationRate);

% Calculate differential rate and median values
differenceRate = avgTargetCreationRate(1:endOfTask) - avgHitRate(1:endOfTask);
medianDifferenceRate = movmedian(differenceRate,1001);
madDifferenceRate = movmad(differenceRate,1001);

%% STEADY STATE RATE
% Determine the last frame for which the lower bound was below 0.1
steadyStateFrame = 1;
for i=endOfTask:-1:1
    if (medianDifferenceRate(i) - madDifferenceRate(i)) < 0.1
        steadyStateFrame = i;
        break
    end
end

% Calculate the steady state rate
steadyTimes = (steadyStateFrame:endOfTask)/200;
overwhelmedSteadyStateRate = mean(hits(steadyStateFrame:endOfTask));

% Catch some edge cases. If the steady state rate is NaN, that's a
% problem. If the steady state rate is 0 then the subject never got
% overwhelmed and we use the maximum hit rate instead.
if isnan(overwhelmedSteadyStateRate)
    minkeyboard
elseif (length(steadyTimes) < 200) || (overwhelmedSteadyStateRate == 0)
    overwhelmedSteadyStateRate = max(avgHitRate);
    steadyStateFrame = endOfTask;
end

% Keep track of all candidate steady state frames for graphing results
candidateSteadyStateFrames = steadyStateFrame;

if verbose
    fprintf('Initial: The Overwhelmed phase lasts from frame %i to %i (%0.2f s to %0.2f s).\n',steadyStateFrame,endOfTask,steadyStateFrame/200,endOfTask/200);
    fprintf('The average hit rate at the start of the phase was %0.2f.\n',avgHitRate(steadyStateFrame));
    fprintf('The median difference rate interval at the start of the phase was (%0.2f, %0.2f).\n',...
        medianDifferenceRate(steadyStateFrame) - madDifferenceRate(steadyStateFrame),medianDifferenceRate(steadyStateFrame) + madDifferenceRate(steadyStateFrame));
    fprintf('The subject achieved a steady state rate of %0.2f.\n',overwhelmedSteadyStateRate);
end

%% OVERWHELMED
% The steadyStateRate is the number of balls contacted
while 1
    % Check to see if overwhelmed phase should start earlier based on the
    % subject performing at a higher level than the calculated steady state
    % rate calculated.
    candidateFrame = steadyStateFrame;
    for i=steadyStateFrame:-1:1
        if (avgHitRate(i) > mean(hits(i:endOfTask)))
            candidateFrame = i;
            candidateSteadyStateRate = mean(hits(i:endOfTask));
            
            % If we are within 5 seconds of starting, give the search a chance
            % chance to keep finding "better" windows.
        elseif i < (steadyStateFrame - 1000)
            break
        end
        
    end
    
    % Record the current candidate frame if it is different
    if candidateSteadyStateFrames(end) ~= candidateFrame
        candidateSteadyStateFrames = [candidateSteadyStateFrames candidateFrame]; %#ok<*AGROW>
    end
    
    if verbose
        fprintf('Higher Performance: The Overwhelmed phase lasts from frame %i to %i (%0.2f s to %0.2f s).\n',candidateFrame,endOfTask,candidateFrame/200,endOfTask/200);
        fprintf('The average hit rate at the start of the phase was %0.2f.\n',avgHitRate(candidateFrame));
        fprintf('The median difference rate interval at the start of the phase was (%0.2f, %0.2f).\n',...
            medianDifferenceRate(candidateFrame) - madDifferenceRate(candidateFrame),medianDifferenceRate(candidateFrame) + madDifferenceRate(candidateFrame));
        fprintf('The subject achieved a steady state rate of %0.2f.\n',candidateSteadyStateRate);
    end
    
    % If the subject was overwhelmed at this new candidate frame, keep
    % pulling back to see where they weren't overhelmed.
    newStartFrame = candidateFrame;
    for i=newStartFrame:-1:1
        if (medianDifferenceRate(i) - madDifferenceRate(i)) > 0.1
            candidateFrame = i;
            candidateSteadyStateRate = mean(hits(candidateFrame:endOfTask));
            
            % If we are within 5 seconds of starting, give the search a chance
            % chance to keep finding "better" windows.
        elseif i < (newStartFrame - 1000)
            break
        end
    end
    
    % Record the current candidate frame if it is different
    if candidateSteadyStateFrames(end) ~= candidateFrame
        candidateSteadyStateFrames = [candidateSteadyStateFrames candidateFrame];
    end
    
    if verbose
        fprintf('Earlier Overwhelmed: The Overwhelmed phase lasts from frame %i to %i (%0.2f s to %0.2f s).\n',candidateFrame,endOfTask,candidateFrame/200,endOfTask/200);
        fprintf('The average hit rate at the start of the phase was %0.2f.\n',avgHitRate(candidateFrame));
        fprintf('The median difference rate interval at the start of the phase was (%0.2f, %0.2f).\n',...
            medianDifferenceRate(candidateFrame) - madDifferenceRate(candidateFrame),medianDifferenceRate(candidateFrame) + madDifferenceRate(candidateFrame));
        fprintf('The subject achieved a steady state rate of %0.2f.\n',candidateSteadyStateRate);
    end
    
    % Update steady state frame if using the candidate frame rewards the
    % subject with a higher steady state rate.
    if candidateFrame < steadyStateFrame
        steadyStateFrame = candidateFrame;
        steadyTimes = (steadyStateFrame:endOfTask)/200;
        overwhelmedSteadyStateRate = candidateSteadyStateRate;
    else
        % Calculate comparison steady state values
        y = quantile(avgHitRate(steadyStateFrame:endOfTask),[0.25 0.50 0.75]);
        overwhelmedQ1 = y(1);
        overwhelmedHitRateMedian = y(2);
        overwhelmedQ3 = y(3);
        overwhelmedHitRateIQR = overwhelmedQ3 - overwhelmedQ1;
        %overwhelmedHitRateMean = mean(avgHitRate(steadyStateFrame:endOfTask));
        steadyStateRate = overwhelmedSteadyStateRate;
        
        if verbose
            fprintf('The subject''s median hit rate during the overwhelmed phase was %0.2f and the hit rate IQR was %0.2f.\n',...
                overwhelmedHitRateMedian, overwhelmedHitRateIQR);
        end
        break
    end
end

%% AT THE LIMIT
% The point at which the subject is at the limit is when the object creation
% rate first reaches their steady state rate or when they stop keeping up
% with the object creation rate leading into the overwhelmed phase, whichever
% point comes first.

% Check for the first frame before the steady state frame at which the
% target creation rate is larger than the subject's steady state rate. Take
% this as the atLimitFrame.
for i=1:steadyStateFrame
    if (avgTargetCreationRate(i) > steadyStateRate)
        break
    end
end
atLimitFrame = i;

% Go back to check whether the subject stops keeping up with the target
% creation rate before the current estimate for the atLimitFrame. If
% yes, take this frame as the new atLimitFrame.
for i=atLimitFrame:-1:1
    if (medianDifferenceRate(i) - madDifferenceRate(i)) > 0.05
        atLimitFrame = i;
    elseif i < steadyStateFrame - 2000
        break
    end
end

% Determine the atLimit interval and the steady state rate within it
atLimitTimes = (atLimitFrame:steadyStateFrame)/200;
if atLimitFrame == steadyStateFrame
    atLimitSteadyStateRate = avgHitRate(atLimitFrame);
else
    atLimitSteadyStateRate = mean(avgHitRate(atLimitFrame:steadyStateFrame));
end

if verbose
    fprintf('The At Limit phase lasts from frame %i to %i.\n',atLimitFrame,steadyStateFrame);
    fprintf('The subject achieved a steady state rate in this phase of %0.2f (trial level steady state rate is %0.2f).\n',atLimitSteadyStateRate,steadyStateRate);
end

% Give subject credit if they have a higher steady state rate in the at
% limit phase than they do in the overwhelmed phase.
if atLimitSteadyStateRate > steadyStateRate
    steadyStateRate = mean(hits(atLimitFrame:endOfTask));
end


%% EASY TIME
% The remaining time in the trial is coded as easy. A subject's baseline
% accuracy is the number of targets hit during this phase.
easyTimes = (1:atLimitFrame)/200;
easyRates = linspace(avgHitRate(1),atLimitSteadyStateRate-0.1,length(easyTimes));
subjectAccuracy = sum(hits(1:atLimitFrame)) / sum(hits(1:atLimitFrame) + misses(1:atLimitFrame));

%% SUMMARIZE RESULTS
if verbose
    fprintf('The phases took place between the following frame numbers.\nEasy: 1-%i, At-Limit: %i-%i, Overwhelmed: %i-%i.\n',...
        atLimitFrame,atLimitFrame,steadyStateFrame,steadyStateFrame,endOfTask);
    fprintf('The subject achieved a steady state rate of %0.2f Hz during the At-Limit phase\nand a steady state rate of %0.2f Hz during the Overwhelmed phase.\n',...
        atLimitSteadyStateRate,overwhelmedSteadyStateRate);
    fprintf('The subject''s overall steady state rate is %0.2f Hz.\n',steadyStateRate);
    fprintf('The subject''s median hit rate for the Overwhelmed phase was %0.2f\nwith an interquartile range of %0.2f.\n',...
        overwhelmedHitRateMedian,overwhelmedHitRateIQR);
end

if graphFlag
    x = (1:endOfTask)/200;
    figure('WindowState','maximized');
    set(gcf,'Visible','on');
    hold on;
    plot(x,avgTargetCreationRate(1:endOfTask),'Color',[0 0 0.8],'LineWidth',4,'DisplayName','Target Creation Rate');
    plot(x,avgHitRate(1:endOfTask),'Color',[0 0.8 0],'LineWidth',4,'DisplayName','Target Contact Rate');
    plot(easyTimes,easyRates,'Color',[0.9290 0.6940 0.1250],'LineWidth',2,'DisplayName','Easy Phase');
    plot(atLimitTimes,atLimitSteadyStateRate*ones(length(atLimitTimes),1),'Color',[0.8500 0.3250 0.0980],'LineWidth',2,'DisplayName','At Limit Phase');
    plot(steadyTimes,overwhelmedSteadyStateRate*ones(length(steadyTimes),1),'Color',[0.6350 0.0780 0.1840],'LineWidth',2,'DisplayName','Overwhelmed Phase');
    
    titleString = sprintf("Steady state rate is %1.2f targets/second",steadyStateRate);
    xl = xlim;
    overwhelmedMidPoint = median(steadyTimes)/xl(2);
    lowHitRate = min(avgTargetCreationRate(steadyStateFrame:endOfTask))/2;
    yl = ylim;
    annotation('textarrow',[overwhelmedMidPoint overwhelmedMidPoint],[lowHitRate/yl(2) steadyStateRate/yl(2)],'String',titleString,'Color',[0.6350 0.0780 0.1840],...
        'LineWidth',2,'FontSize',24);
    lgd = legend('location','northwest');
    xlabel('Time (s)');
    ylabel('Rates (Hz)');
    grid on;
    ax = gca;
    ax.FontSize = 32;
    ax.FontWeight = 'bold';
    ax.LineWidth = 2.5;
    lgd.FontSize = 24;
    lgd.FontWeight = 'normal';
    hold off;
end

if steadyStateRate > 5
    keyboard
end

end
