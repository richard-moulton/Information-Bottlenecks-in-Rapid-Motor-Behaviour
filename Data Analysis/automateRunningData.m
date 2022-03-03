%% A script to automate running all OH, OHA, TOH, and TOHA trials
%  The "tree" of function calls is for this script to call
%  processMultipleOHFiles.m for each block of trials in a task; this calls
%  processSingleOHFile.m for each trial in the block; and this calls the
%  necessary series of function to convert a continuous trial into an
%  event-based representation and analyse it.
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
%  3 December, 2021

% Blank slate start
clear variables;
close all;
clc;

% File-level settings
graphFlag = 0;      % If 1, produce figures throughout
diaryFlag = 0;      % If 1, print decision list data out to file
verbose = 0;        % If 1, print everything to screen.
filesFlag = 0;      % 3 full list (OH), 5 full list (OHA), 6-7 turbo variants of tasks, 10-11 are OH and TOH comparisons
maxNumRecords = 500;% Maximum size for a block

addpath('../Utilities/');
fileSaveType = 'eps';% Format in which to save figures

% Record the names of the files produced
ohListNames = {};
ohaListNames = {};
tohListNames = {};
tohaListNames = {};

% Analyse all Object-Hit trials
startDir = 1;       % Directory number at which to start
while 1
    [endFlag, endDir, matFileName] = processMultipleOHFiles(graphFlag,diaryFlag,3,verbose,startDir,maxNumRecords,fileSaveType);
    
    ohListNames{length(ohListNames)+1} = matFileName;
    
    if endFlag
        break
    else
        startDir = endDir +1;
    end    
end

% Follow up by analysing those subjects who completed OH and TOH
% There are only 38 subjects, so no combining is required
startDir = 1;       % Directory number at which to start
processMultipleOHFiles(graphFlag,diaryFlag,10,verbose,startDir,maxNumRecords,fileSaveType);

% Analyse all Object-Hit-and-Avoid trials
startDir = 1;       % Directory number at which to start
while 1
    [endFlag, endDir, matFileName] = processMultipleOHFiles(graphFlag,diaryFlag,5,verbose,startDir,maxNumRecords,fileSaveType);
    
    ohaListNames{length(ohaListNames)+1} = matFileName;
    
    if endFlag
        break
    else
        startDir = endDir + 1;
    end    
end

% Analyse all Turbo Object-Hit trials
startDir = 1;       % Directory number at which to start
while 1
    [endFlag, endDir, matFileName] = processMultipleOHFiles(graphFlag,diaryFlag,6,verbose,startDir,maxNumRecords,fileSaveType);
    
    tohListNames{length(tohListNames)+1} = matFileName; %#ok<*SAGROW>
    
    if endFlag
        break
    else
        startDir = endDir + 1;
    end    
end

% Follow up by analysing those subjects who completed TOH and OH
% There are only 38 subjects, so no combining is required
startDir = 1;       % Directory number at which to start
processMultipleOHFiles(graphFlag,diaryFlag,11,verbose,startDir,maxNumRecords,fileSaveType);

% Turbo Object-Hit-and-Avoid
startDir = 1;       % Directory number at which to start
while 1
    [endFlag, endDir, matFileName] = processMultipleOHFiles(graphFlag,diaryFlag,7,verbose,startDir,maxNumRecords,fileSaveType);
    
    tohaListNames{length(tohaListNames)+1} = matFileName;
    
    if endFlag
        break
    else
        startDir = endDir + 1;
    end    
end

matFileCombiner(ohListNames,ohaListNames,tohListNames,tohaListNames);
