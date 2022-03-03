%% Extract state-space representations for each event in an OH/TOH/OHA/TOHA trial
%  The generated state-space representation includes paddle positions as
%  well as object positions, speed and values. Object IDs are included to 
%  make trial-level tracking of objects easier.
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

function statesList = getStatesList(eventsList,objectsList,data)

statesList = cell(size(eventsList,1),1);        % an empty array for storing the states that result from the interpreted events.

% For each event that occurs, calculate the state that it leads to.
for ii=1:length(eventsList)
    statesList(ii) = {getState(eventsList(ii,:),objectsList,data)};
end

end

%% Functions required to run getStatesList

% Generate and return the state vector for a given frame. Note, while this
% could be called for every single frame, the idea is to only compute these
% when events occur.
function [state] = getState(event,objectsList,data)
paddleStates = getPaddleStates(data,event(1));
objStates = getObjStates(data,event,objectsList);

state = cell([paddleStates objStates]);
end

% Generate the portion of the state vector relating to the position of the
% paddles in the workspace.
function paddleStates = getPaddleStates(data,frame)
paddle1X = data.c3d.Left_HandX(frame,1);
paddle1Y = 100*data.c3d.Left_HandY(frame,1);
paddle2X = data.c3d.Right_HandX(frame,1);
paddle2Y = 100*data.c3d.Right_HandY(frame,1);

paddleStates = {[paddle1X paddle1Y] [paddle2X paddle2Y]};
end

% Generate and return the portion of the state vector that relates to the
% position of objects in the workspace.
function [objStates] = getObjStates(data,event,objectsList)
frame = event(1);
numObjs = 1;
objStates = cell(1);

for ii = 1:length(objectsList)
    if (objectsList(ii).getEndFrame < frame)
        continue
    elseif(event(2) == 2) && (event(3) == objectsList(ii).getID())
        continue
    elseif(objectsList(ii).getStartFrame > frame)
        break
    elseif ((objectsList(ii).getStartFrame() <= frame) && (objectsList(ii).getEndFrame() > frame))
        [~,yCH] = getXYbyChannel(data,objectsList(ii).getChannel());
        objectY = yCH(frame);
        
        currObject = ohObject(objectsList(ii).getID(),objectsList(ii).getBin(),100*objectY,objectsList(ii).getSpeed(),objectsList(ii).getValue());
        
        objStates{numObjs} = currObject;
        numObjs = numObjs + 1;
    end
end
end
