%% An object (target or distractor) in OH/OHA/TOH/TOHA
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
%  This class stores trial-level information about an object in a single
%  OH/OHA/TOH/TOHA trial.
%
%  For a snapshot in time of a single object, see the ohObject class.
%
%  6 December, 2021

classdef trialObject
    
    properties
        % These properties have the nonsensical NaN values as their default
        % value. This can be used to see whether the object has actual data
        % associated with it or if it is simply a preallocation of memory.
        id = NaN                      % the object's position in the sequential order
        value = NaN                   % the gain/loss from hitting the object
        channel = NaN                 % the c3d channel to which this object is assigned
        bin = NaN                     % the bin to which this object is assigned
        speed = NaN                   % the object's initial speed, calculated from its initial trajectories
        startFrame = NaN              % the frame in which this object appeared
        contactFrame = NaN            % the frame in which this object was hit
        endFrame = NaN                % the frame in which this object disappears
        endReason = NaN               % 1 if the object was successfully hit, larger if it was missed (different kinds of misses are encoded differently)
                
        % These arrays contain the object's complete trajectory
        xTrajectory = [];
        yTrajectory = [];
    end
    
    methods
        % CONSTRUCTOR numObjects,distractorFlag,channel,bin,startFrame
        function obj = trialObject(id,val,ch,bN,sF)
            if nargin > 0
                obj.id = id;
                obj.value = val;
                obj.channel = ch;
                obj.bin = bN;
                obj.startFrame = sF;
            end
        end
        
        % GETTERS
        function val = getID(obj)
            val = obj.id;
        end
        function val = getValue(obj)
           val = obj.value; 
        end
        function val = getChannel(obj)
            val = obj.channel;
        end
        function val = getBin(obj)
            val = obj.bin;
        end
        function val = getSpeed(obj)
            val = obj.speed;
        end
        function val = getStartFrame(obj)
            val = obj.startFrame;
        end
        function val = getContactFrame(obj)
           val = obj.contactFrame; 
        end
        function val = getEndFrame(obj)
            val = obj.endFrame;
        end
        function val = getEndReason(obj)
            val = obj.endReason;
        end
        function hitFlag = wasHit(obj)
            if isnan(obj.contactFrame)
                hitFlag = false;
            else
                hitFlag = true;
            end
        end
        
        % SETTERS
        function obj = setBin(obj,bin)
            obj.bin = bin;
        end
        function obj = setContactFrame(obj,frame)
            obj.contactFrame = frame;
        end
        function obj = setTrajectories(obj,xTraj,yTraj,bounds)
            % Assign the trajectory arrays
            obj.xTrajectory = xTraj;
            obj.yTrajectory = yTraj;
            
            % Calculate the object's initial speed in terms of y-bins/s
            binHeight = (bounds.yMaximum - bounds.yMinimum) / 10;
            objectSpeed = (yTraj(2) - yTraj(1)) * 2;
            obj.speed = objectSpeed / binHeight;         
            
            % Note the object's end frame
            obj.endFrame = obj.startFrame + length(obj.xTrajectory) - 1;
            
            % Reason about the object's end
            if obj.wasHit()
                if obj.yTrajectory(end) >= obj.yTrajectory(end-1)
                    obj.endReason = 1;       % If the object left the screen and wasn't going down, then it was hit
                elseif obj.xTrajectory(end) == obj.xTrajectory(end-1)
                    obj.endReason = 2;       % If the object is exiting in a straight line then it was missed.
                else
                    
                    % If the object left the screen while moving down on an
                    % angle, we have to recover it's final path to judge
                    % whether or not it exited out the bottom of the screen.
                    slope = (obj.yTrajectory(end) - obj.yTrajectory(end-1)) / (obj.xTrajectory(end) - obj.xTrajectory(end-1));
                    intercept = obj.yTrajectory(end) - (slope * obj.xTrajectory(end));
                    
                    if slope > 0
                        if  slope * bounds.getXMin() + intercept <= bounds.getYMin()
                            obj.endReason = 3;   % If the path didn't get to the edge of the screen before dropping off, the object was missed
                        else
                            obj.endReason = 1;
                        end
                    else
                        if  slope * bounds.getXMax() + intercept <= bounds.getYMin()
                            obj.endReason = 3;   % If the path didn't get to the edge of the screen before dropping off, the object was missed
                        else
                            obj.endReason = 1;
                        end
                    end
                end
            else
                obj.endReason = 4;       % If the object was never hit, then it was missed
            end
            
        end
    end
end
