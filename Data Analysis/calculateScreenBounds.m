%% Calculate x- and y-bounds for an OH/OHA/TOH/TOHA trial's workspace
%  Determine the minimum and maximum values in the x- and y-dimensions of
%  the workspace for an OH/OHA/TOH/TOHA task.
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

function bounds = calculateScreenBounds(data)
xCentre = (data.c3d.CALIBRATION.LEFT_SHO_X + data.c3d.CALIBRATION.RIGHT_SHO_X)/2;
trialProtocol = data.c3d.TRIAL.TP;
screenFrame = data.c3d.TP_TABLE.Screen_Frame(trialProtocol);
yCentre = data.c3d.TARGET_TABLE.Y_GLOBAL(screenFrame);
screenHeight = data.c3d.TARGET_TABLE.Height(trialProtocol);
screenWidth = data.c3d.TARGET_TABLE.Width(trialProtocol);

xMin = xCentre - (screenWidth/2);
xMax = xCentre + (screenWidth/2);
yMin = yCentre - (screenHeight/2);
yMax = yCentre + (screenHeight/2);

bounds = screenBounds(xMin,xMax,yMin,yMax);
end
