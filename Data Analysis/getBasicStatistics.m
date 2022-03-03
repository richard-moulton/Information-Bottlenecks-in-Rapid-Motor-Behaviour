%% Calculate basic statistics for an OH/OHA/TOH/TOHA trial
%  For a single trial, calculate the frames where targets and distractors
%  are created and the frames where hits, misses, drops (successfully
%  missed distractors), and crashes (unsuccessfully avoided distractors)
%  occur.
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

function [targetCreates,distractorCreates,hits,misses,drops,crashes] = getBasicStatistics(data,eventsList)

numberOfFrames = length(data.c3d.event_code);

% Initialize arrays for basic statistics
targetCreates = zeros(numberOfFrames,1);
distractorCreates = zeros(numberOfFrames,1);
hits = zeros(numberOfFrames,1);
misses = zeros(numberOfFrames,1);
drops = zeros(numberOfFrames,1);
crashes = zeros(numberOfFrames,1);

% Trial-level summary statistics
numHits = 0;
numMisses = 0;
numDrops = 0;
numCrashes = 0;
numEvents = 1;

% This array is used to track objects' fates
idTracker = zeros(300,1);

% Perform analysis on the basis of frames. This include counting hits and
% misses as the subject progresses through the trial, tracking object
% creation and contact rates, the centre-of-gravity of the balls on the
% screen and the total number of balls being faced.
for frameNum = 1:numberOfFrames
    
    while eventsList(numEvents,1) == frameNum
        switch eventsList(numEvents,2)
            case 1      % Ball Creation
                targetCreates(frameNum) = 1;
            case 2
            case 3      % Ball End
                switch eventsList(numEvents,4)
                    case 1
                        numHits = numHits + 1;
                        hits(frameNum) = hits(frameNum) + 1;
                        idTracker(eventsList(numEvents,3)) = idTracker(eventsList(numEvents,3)) + 1;
                    otherwise
                        numMisses = numMisses + 1;
                        misses(frameNum) = misses(frameNum) + 1;
                        idTracker(eventsList(numEvents,3)) = idTracker(eventsList(numEvents,3)) + 1;
                end
            case 4
                distractorCreates(frameNum) = misses(frameNum) + 1;
            case 5
                numCrashes = numCrashes + 1;
                crashes(frameNum) = crashes(frameNum) + 1;
                idTracker(eventsList(numEvents,3)) = idTracker(eventsList(numEvents,3)) + 1;
            case 6
                switch eventsList(numEvents,4)
                    case 4
                        numDrops = numDrops + 1;
                        drops(frameNum) = drops(frameNum) + 1;
                        idTracker(eventsList(numEvents,3)) = idTracker(eventsList(numEvents,3)) + 1;
                end
            otherwise
                fprintf("Came across a poorly formatted event: [%i %i %i %i %i].\n",eventsList(numEvents,1),eventsList(numEvents,2),eventsList(numEvents,3),eventsList(numEvents,4),eventsList(numEvents,5));
        end
        
        % Done processing this event, get ready for the next
        numEvents = numEvents + 1;
        
        if numEvents > length(eventsList)
            break
        end
    end
    
    if numEvents > length(eventsList)
        break
    end
end

% Trim the returned arrays if necessary
if frameNum < numberOfFrames
    targetCreates = targetCreates(1:frameNum);
    distractorCreates = distractorCreates(1:frameNum);
    hits = hits(1:frameNum);
    misses = misses(1:frameNum);
    drops = drops(1:frameNum);
    crashes = crashes(1:frameNum);
end

% Error-checking.
if numHits + numMisses + numDrops + numCrashes ~= 300
    keyboard
end
end
