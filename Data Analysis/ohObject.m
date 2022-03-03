%% An object in OH/OHA/TOH/TOHA at a specific moment in time.
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
%  This class stores all of the recorded information about a single object
%  at a single, specific frame during a trial.
%
%  For trial-level information about an object, see the trialObject class.
%
%  6 December, 2021

classdef ohObject
    
    properties
        id = -1;    % The object's ID (useful for cross-referencing representations)
        xBin = -1;  % The bin for the object's location in the x-axis
        yBin = -1;  % The bin for the object's location in the y-axis
        speed = -1; % The object's speed
        value = 0; % A flag indicating the object's value.
                    % (OH: all objects are worth 1. OHA: objects 1,
                    % distractors -1).
    end
    
    methods
        %CONSTRUCTOR
        function obj = ohObject(id,x,y,speed,value)
           obj.id = id;
           obj.xBin = x;
           obj.yBin = y;
           obj.speed = speed;
           obj.value = value;
        end
        
        % GETTERS
        function val = getID(obj)
            val = obj.id;
        end
        function val = getXBin(obj)
            val = obj.xBin;
        end
        function val = getYBin(obj)
            val = obj.yBin;
        end
        function val = getSpeed(obj)
            val = obj.speed;
        end
        function val = getValue(obj)
            val = obj.value;
        end
    end
end
