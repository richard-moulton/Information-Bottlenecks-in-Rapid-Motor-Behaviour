%% Determine when decisions are made during an OH/OHA/TOH/TOHA trial
%  Produces lists of decisions in an OH/OHA/TOH/TOHA trial based on two
%  criteria: actions represent decisions, so any contact with an object is
%  evidence that a decision was made; and that decisions are made on state
%  changes in the trial.
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

function [decisionList,decisionFrames,gaps] = getDecisionList(eventsList,statesList,atLimitFrame)

% Isolate the frame number data from eventsList
frames = eventsList(:,1);

% Initialize memory
decisionList = cell(size(eventsList,1),1);        % a (too large) empty array for storing decisions
decisionFrames = NaN(size(eventsList,1),1);      % a (too large) empty array for recording at which frames decision are made
numDecisions = 0;

% Sanity check that we can find the desired start event
if(atLimitFrame < 1) || (atLimitFrame > frames(end))
    fprintf("ALERT!! Looks like we have an improperly calculated frame number: %i (outside the range [%i, %i]).\n",atLimitFrame,1,frames(end));
    keyboard
end

startEventIndex = findStartEventIndex(frames,1,length(frames),atLimitFrame);
gaps = zeros(size(eventsList,1),1);

for eventNum = startEventIndex:length(eventsList)
    if eventsList(eventNum,2) == 2
        ballID = eventsList(eventNum,3);
        
        paddleNum = eventsList(eventNum,4);
        [decisionState,decisionFrame] = getDecisionState(eventNum,statesList,eventsList,ballID);
        
        numDecisions = numDecisions+1;
        decisionList(numDecisions) = {[paddleNum ballID decisionState]};
        decisionFrames(numDecisions) = decisionFrame;
        gaps(numDecisions) = eventsList(eventNum,1) - decisionFrame;
        
        if gaps(numDecisions) < 0
            fprintf("WARNING!! The gap is somehow negative for event %i!\n",eventNum);
            keyboard
        end
    end
end

decisionList = decisionList(1:numDecisions);
decisionFrames = decisionFrames(1:numDecisions);
gaps = gaps(1:numDecisions);

end

%% Functions required to run getDecisionList

% Implement a binary search to find the index for the first event after the
% argument frame number (used here to find an event index relative to the
% atLimitFrame).
function startIndex = findStartEventIndex(frames,first,last,startFrame)

if first == last
    startIndex = first;
elseif frames(1) > startFrame
    startIndex = 1;
else
    midpoint = floor(first + (last - first)/2);
    
    if midpoint == first
        startIndex = last;
    elseif(frames(midpoint) < startFrame)
        startIndex = findStartEventIndex(frames,midpoint,last,startFrame);
    elseif(frames(midpoint) >= startFrame)
        startIndex = findStartEventIndex(frames,first,midpoint,startFrame);
    else
        fprintf("ALERT!! It is possible that multiple identical frame numbers are interfering with the binary search!");
    end
end
end

% Find the state that happened after the most recent paddle contact (and at
% least before 100 ms (20 frames) before the argument frame number). This
% represents the necessary lag between visual stimulus and voluntary action
% [CITATION NEEDED] and potentially lines up better with a sequential
% decision making model.
function [decisionState,decisionFrame] = getDecisionState(startEvent,statesList,eventsList,ballID)
startFramePoint = eventsList(startEvent,1);
targetFramePoint = startFramePoint - 20;

% Working backwards from the start event, we look for a ball contact event
% that occured before the target frame point.
for ii=startEvent:-1:1
    if (eventsList(ii,1) < targetFramePoint) && (eventsList(ii,2) == 2)
        if isInState(statesList{ii},ballID)
            decisionState = statesList{ii};
            decisionFrame = eventsList(ii,1);
            return
        else
            break
        end
    end
end

% If we can't, then working forwards we look for the event where the
% contacted ball was created.
for jj=ii:startEvent-1
    if (eventsList(jj,2) == 1) && (isInState(statesList{jj},ballID))
        decisionState = statesList{jj};
        decisionFrame = eventsList(jj,1);
        return
    end
end

end

% Returns true if the argument ballID is present in the argument state
function flag = isInState(state,ballID)
flag = false;
objects = state(3:end);

for i=1:length(objects)
    obj = objects(i);
    if isempty(obj{1})
        continue
    elseif(obj{1}.getID() == ballID)
        flag = true;
        return
    end
end
end
