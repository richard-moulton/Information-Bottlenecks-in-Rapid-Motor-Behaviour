%% This function calculates a subject's accuracy on a per-bin basis
%  Subject accurady is divided between the easy/at-limit phases and the
%  overwhelmed phase. Targets and distractors are considered separately.
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

function [earlyBinAccuracy,overwhelmedBinAccuracy,earlyBinAvoidance,overwhelmedBinAvoidance] = getBinAccuracy(eventsList,objectsList,steadyStateFrame,verbose)
% Initialize arrays
earlyBinHits = zeros(10,1);
earlyBinTargets = zeros(10,1);
earlyBinMisses = zeros(10,1);
earlyBinDistractors = zeros(10,1);

overwhelmedBinHits = zeros(10,1);
overwhelmedBinObjects = zeros(10,1);
overwhelmedBinMisses = zeros(10,1);
overwhelmedBinDistractors = zeros(10,1);

xBins = unique(eventsList(eventsList(:,2)==1,4));
xBinsMap = containers.Map(xBins,linspace(1,length(xBins),length(xBins)));

% Go through each object and assess whether it was successfully hit or not.
for i=1:length(objectsList)
    bin = xBinsMap(objectsList(i).getBin());
    if objectsList(i).getValue() > 0 % Targets
        if objectsList(i).getStartFrame() < steadyStateFrame
            earlyBinTargets(bin) = earlyBinTargets(bin) + 1;
            if objectsList(i).getEndReason() == 1
                earlyBinHits(bin) = earlyBinHits(bin) + 1;
                if verbose
                    fprintf('Object %i (Bin %i): Target. Hit.\n',i,bin);
                end
            else
                if verbose
                    fprintf('Object %i (Bin %i): Target. Missed.\n',i,bin);
                end
            end
        else
            overwhelmedBinObjects(bin) = overwhelmedBinObjects(bin) + 1;
            if objectsList(i).getEndReason() == 1
                overwhelmedBinHits(bin) = overwhelmedBinHits(bin) + 1;
            end
        end
    else % Distractors
        if objectsList(i).getStartFrame() < steadyStateFrame
            earlyBinDistractors(bin) = earlyBinDistractors(bin) + 1;
            if objectsList(i).getEndReason() == 4
                earlyBinMisses(bin) = earlyBinMisses(bin) + 1;
                if verbose
                    fprintf('Object %i (Bin %i): Distractor. Dropped.\n',i,bin);
                end
            else
                if verbose
                    fprintf('Object %i (Bin %i): Distractor. Crashed.\n',i,bin);
                end
            end
        else
            overwhelmedBinDistractors(bin) = overwhelmedBinDistractors(bin) + 1;
            if objectsList(i).getEndReason() == 4
                overwhelmedBinMisses(bin) = overwhelmedBinMisses(bin) + 1;
            end
        end
    end
end

% Calculate accuracies
earlyBinAccuracy = earlyBinHits./earlyBinTargets;
overwhelmedBinAccuracy = overwhelmedBinHits./overwhelmedBinObjects;

earlyBinAvoidance = earlyBinMisses./earlyBinDistractors;
overwhelmedBinAvoidance = overwhelmedBinMisses./overwhelmedBinDistractors;
end
