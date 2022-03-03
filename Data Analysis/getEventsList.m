%% Convert "continuous" data from OH/OHA/TOH/TOHA into events and objects.
%  By going through the frames that have been captured at 200 Hz, this
%  function extracts the relevant events and objects to reduce file size
%  while allowing all the necessary analysis.
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

function [eventsList,objectsList] = getEventsList(data,NUM_OBJS)
% Initialize arrays and variables
bounds = calculateScreenBounds(data);       % record the screen's bounds based on the calibration data that has been loaded into 'data'
eventIdx = find(data.c3d.event_code);       % indices for events throughout the trial. Useful for addressing into any other array for synchronized readings.
eventsList = [];                            % an empty array for storing interpreted events.
numEvents = 0;                              % a counter for the number of events seen to date.
objectsList(NUM_OBJS) = trialObject;        % a list of the objects from the task
numObjects = 0;                             % a counter for the number of objects seen to date

% Rename the event-related data structures that will be regularly referenced.
eventCodes = abs(data.c3d.event_code);
eventSupplements = data.c3d.event_supplement;

% Extract and interpret all events during a trial
for ii = 1:length(eventIdx)
    % An event code of 100 or -100 indicates the creation of an object
    if(eventCodes(eventIdx(ii)) == 100)
        [value, channel, ~, ~] = getObjCreationSupplementalDetails(eventSupplements(eventIdx(ii)));
        startFrame = eventIdx(ii);
        
        [ballCH_x,~] = getXYbyChannel(data,channel);
        bin = ballCH_x(startFrame);
        
        % Create a trialObject to track
        numObjects = numObjects + 1;
        objectsList(numObjects) = trialObject(numObjects,value,channel,bin,startFrame);
        
        % The new event is coded by the frame in which it occurred, its
        % number and name, a flag indicating if the object is a distractor
        % or not, the channel the object was assigned to, the object's
        % targetTableRow (?) and the remainder of the bit code (should
        % be zeros (?))
        numEvents = numEvents + 1;
        eventsList = [eventsList;[startFrame, 1, numObjects, bin, value]]; %#ok<*AGROW>
        
        % All other event codes indicate a paddle contacting a ball
    else
        [contactPaddle1, contactPaddle1b, contactPaddle2, contactPaddle2b] = getPaddleContactSupplementalDetails(eventSupplements(eventIdx(ii)));
        frameNum = eventIdx(ii);

        % Determine which paddle made contact with the object
        if (contactPaddle1 > 0)
            contactPaddle = 1;
            
            [newEvent, objectsList] = getObjContactEvent(frameNum,objectsList,contactPaddle,contactPaddle1);
            numEvents = numEvents + 1;
            eventsList = [eventsList;newEvent];
        end
        if (contactPaddle1b > 0)
            contactPaddle = 1;
            
            [newEvent, objectsList] = getObjContactEvent(frameNum,objectsList,contactPaddle,contactPaddle1b);
            numEvents = numEvents + 1;
            eventsList = [eventsList;newEvent];
        end
        if (contactPaddle2 > 0)
            contactPaddle = 2;
            
            [newEvent, objectsList] = getObjContactEvent(frameNum,objectsList,contactPaddle,contactPaddle2);
            numEvents = numEvents + 1;
            eventsList = [eventsList;newEvent];
        end
        if (contactPaddle2b > 0)
            contactPaddle = 2;
            
            [newEvent, objectsList] = getObjContactEvent(frameNum,objectsList,contactPaddle,contactPaddle2b);
            numEvents = numEvents + 1;
            eventsList = [eventsList;newEvent];
        end
        if ((contactPaddle1 == 0) && (contactPaddle1b == 0) && (contactPaddle2 == 0) && (contactPaddle2b == 0))
            fprintf("ALERT! We have an event supplement encoding contact with neither paddle! Frame %f, [%f %f %f %f].\n",frameNum,contactPaddle1, contactPaddle1b, contactPaddle2, contactPaddle2b);
        end
    end
end

% Iterate through all balls in object list to assign each ball its
% trajectory and to generate events representing the ball's ultimate fate.
for ii=1:length(objectsList)
    if isnan(objectsList(ii).getBin())
        fprintf("\nWhoops, a NaN from %s\n",data.file_name);
        fprintf("Waiting...\n");
    end
    
    % Assign trajectories
    objectsList(ii) = setTrajectories(objectsList(ii),data,bounds);
    
    % Create events reflecting the balls' ends
    eventsList= [eventsList;[objectsList(ii).getEndFrame(), 3, objectsList(ii).getID(),objectsList(ii).getEndReason(),objectsList(ii).getValue()]];
end

% If the file label contains "avoid" then it is an OHA task and there are 
% distractors present. Adjust the relevant events.
if contains(data.file_label,["avoid","Avoid"])
    for ii=1:length(eventsList)
       if (objectsList(eventsList(ii,3)).getValue() == -1)
          eventsList(ii,2) = eventsList(ii,2) + 3; 
       end
    end
end

eventsList = sortrows(eventsList);

end

%% Functions required for getEventsList to run.

% Extract variables from the bitwise representation of the event_supplement
% to Obj Creation events.
function [value, channel, targetTableRow, remainder] = getObjCreationSupplementalDetails(eventSupplement)
objectFlag = bi2de(bitget(eventSupplement,1));
if objectFlag
    value = 1;
else
    value = -1;
end
channel = bi2de(bitget(eventSupplement,2:5));
targetTableRow = bi2de(bitget(eventSupplement,6:13));
remainder = bi2de(bitget(eventSupplement,14:16));
end

% Extract variables from the bitwise representation of the event_supplement
% to Paddle Contact events.
function [contactPaddle1, contactPaddle1b, contactPaddle2, contactPaddle2b] = getPaddleContactSupplementalDetails(eventSupplement)
contactPaddle1 = bi2de(bitget(eventSupplement,9:12));
contactPaddle1b = bi2de(bitget(eventSupplement,13:16));
contactPaddle2 = bi2de(bitget(eventSupplement,1:4));
contactPaddle2b = bi2de(bitget(eventSupplement,5:8));
end

% The new event is coded by the frame in which it occurred, the
% paddle involved and name, and then the raw supplemental
% codes. A non-zero supplemental code indicates the channel of
% the contact object
function [newEvent,objectList] = getObjContactEvent(frameNum,objectList,contactPaddle,channel)

% Update the appropriate ball's contact frame. This will ensure
% that each ball stores it's last contact with either paddle.
objContacted = getObjContacted(objectList,frameNum,channel);

% Only create a new event if the ball hasn't been contacted before. This
% will ensure that double hits are not included in the analysis.
if objectList(objContacted).wasHit()
    newEvent = [];
else
    newEvent = [frameNum,2,objContacted,contactPaddle,objectList(objContacted).getValue()];
    objectList = updateContactFrame(objectList,objContacted,frameNum);
end
end

% Based on the argument frame and channel, determine the ID of the object
% that was contacted in that frame.
function objID = getObjContacted(objectList,frame,channel)
objID = -1;
for ii = length(objectList):-1:1
    currObject = objectList(ii);
    if((currObject.getChannel() == channel) && (currObject.getStartFrame() < frame))
        objID = currObject.getID();
        break
    end
end
end

% Update the argument object's contact frame.
function objectList = updateContactFrame(objectList,objID,frame)
currObject = objectList(objID);
objectList(objID) = currObject.setContactFrame(frame);
end

% Assigns x- and y-trajectories to the argument object.
function obj = setTrajectories(object,data,bounds)

% We need to retrieve the object's x- and y-channels
[xCH,yCH] = getXYbyChannel(data,object.getChannel());

% We need to determine which frame in the channel is the object's last.
startFrame = object.getStartFrame();
endFrame = object.getStartFrame();
for ii = object.getStartFrame():length(xCH)
    if xCH(ii) == -1000
        endFrame = ii - 1;
        break
    elseif ii == length(xCH)
        endFrame = ii;
    end
end

% Now we can assign the object its trajectories from the correct portion of
% the correct channel readout.
obj = object.setTrajectories(100*xCH(startFrame:endFrame),100*yCH(startFrame:endFrame),bounds);
end
