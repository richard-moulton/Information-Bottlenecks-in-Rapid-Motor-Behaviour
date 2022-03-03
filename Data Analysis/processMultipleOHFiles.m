%% A function for processing multiple OH/OHA/TOH/TOHA trials.
%  If varargin is not empty, then it must be the
%  variables graphFlag,diaryFlag,filesFlag,verbose,startDir,maxNumRecords.
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
%  Requires: Image Processing Toolbox
%            Statistics and Machine Learning Toolbox
%            Communications Toolbox

function [endFlag, endDir, matFileName] = processMultipleOHFiles(varargin)

addpath('../Utilities/');

if isempty(varargin)
    graphFlag = 0;      % If 1, produce figures throughout
    diaryFlag = 0;      % If 1, print decision list data out to file
    filesFlag = 4;      % 1 single file, 2 small-scale trial (OH), 3 full list (OH), 4 small-scale trial (OHA), 5 full list (OHA), 6-7-8-9 turbo variants of tasks, 10-11 are OH and TOH comparisons
    verbose = 0;        % If 1, print everything to screen.
    
    startDir = 1;       % Directory number at which to start
    maxNumRecords = 500;% Maximum size for a block
    
    fileSaveType = 'png';% Format in which to save figures
else
    graphFlag = varargin{1};
    diaryFlag = varargin{2};
    filesFlag = varargin{3};
    verbose = varargin{4};
    
    startDir = varargin{5};
    maxNumRecords = varargin{6};
    
    fileSaveType = varargin{7};
end

% Where to put trouble records; a) for easy inspection, and b) to avoid
% processing them over and over again if they are corrupt.
troubleRecordsFolder = '/home/richard/Documents/School/PhD/Scott Lab/KINARM data/Trouble Records';

% The filesFlag tells the function which folder to iterate through.
switch filesFlag
    case 2
        OHflag = 1;
        rootDirectory = '/home/richard/Documents/School/PhD/Scott Lab/KINARM data/Small-Scale Data for Trials/OH/';
        matFileName = 'decisionListsOH.mat';
        turboFlag = 0;
    case 3
        OHflag = 1;
        rootDirectory = '/home/richard/Documents/School/PhD/Scott Lab/KINARM data/Classic OH-OHA (control)/OH';
        matFileName = 'fullDecisionListsOH.mat';
        turboFlag = 0;
    case 4
        OHflag = 0;
        rootDirectory = '/home/richard/Documents/School/PhD/Scott Lab/KINARM data/Small-Scale Data for Trials/OHA/';
        matFileName = 'decisionListsOHA.mat';
        turboFlag = 0;
    case 5
        OHflag = 0;
        rootDirectory = '/home/richard/Documents/School/PhD/Scott Lab/KINARM data/Classic OH-OHA (control)/OHA';
        matFileName = 'fullDecisionListsOHA.mat';
        turboFlag = 0;
    case 6
        OHflag = 1;
        rootDirectory = '/home/richard/Documents/School/PhD/Scott Lab/KINARM data/Turbo OH-OHA (Winsport)';
        matFileName = 'fullDecisionListsTurboOH.mat';
        turboFlag = 1;
    case 7
        OHflag = 0;
        rootDirectory = '/home/richard/Documents/School/PhD/Scott Lab/KINARM data/Turbo OH-OHA (Winsport)';
        matFileName = 'fullDecisionListsTurboOHA.mat';
        turboFlag = 1;
    case 8
        OHflag = 1;
        rootDirectory = '/home/richard/Documents/School/PhD/Scott Lab/KINARM data/Small-Scale Data for Trials/Winsport/';
        matFileName = 'decisionListsTurboOH.mat';
        turboFlag = 1;
    case 9
        OHflag = 0;
        rootDirectory = '/home/richard/Documents/School/PhD/Scott Lab/KINARM data/Small-Scale Data for Trials/Winsport/';
        matFileName = 'decisionListsTurboOHA.mat';
        turboFlag = 1;
    case 10
        OHflag = 1;
        rootDirectory = '/home/richard/Documents/School/PhD/Scott Lab/KINARM data/Classic-Turbo OH';
        matFileName = 'decisionListsCompOH.mat';
        turboFlag = 0;
    case 11
        OHflag = 1;
        rootDirectory = '/home/richard/Documents/School/PhD/Scott Lab/KINARM data/Classic-Turbo OH';
        matFileName = 'decisionListsCompTOH.mat';
        turboFlag = 1;
end

fileToLoad = '*.zip';
subjectDirectories = dir2(rootDirectory);
filesSeen = 0;
dirFilesSeen = 0;
totalSubjects = 0;

% Break up large numbers of records into manageable sizes for memory
if startDir ~= 1
    blockFlag = true;
else
    blockFlag = false;
end


% Initialize a structure array for the motor control meta data from each
% trial. Source is data.analysis.REPORT_FEATURES from the c3d file.
if OHflag
    metaDataArray = setObjectHitMetaData();
    taskName = 'Object-Hit';
else
    metaDataArray = setObjectHitAndAvoidMetaData();
    taskName = 'Object-Hit-and-Avoid';
end

% Add the Turbo prefix if necessary.
if turboFlag
    taskName = strcat('Turbo-',taskName);
end

% Estimate an appropriate size for the arrays to store.
memorySize = ceil(1.1* maxNumRecords);

% Initialize arrays
subjectNames = strings(memorySize,1);
trialNames = strings(memorySize,1);
decisionLists = cell(memorySize,1);
decisionFramesList = cell(memorySize,1);
metaDataArray(memorySize).DISTRACTOR_HITS_LEFT = -1;
steadyStateRates = NaN(memorySize,1);
owhrMedian = NaN(memorySize,1);
owhrIQR = NaN(memorySize,1);
steadyStateFrames = NaN(memorySize,1);
atLimitFrames = NaN(memorySize,1);
hitsScores = NaN(memorySize,1);
dropsScores = NaN(memorySize,1);
avgGaps = NaN(memorySize,1);
maxAvgTargetCreationRate = NaN(memorySize,1);
subjectAccuracies = NaN(memorySize,1);
earlyBinAccuracies = NaN(memorySize,10);
overwhelmedBinAccuracies = NaN(memorySize,10);
earlyBinAvoidances = NaN(memorySize,10);
overwhelmedBinAvoidances = NaN(memorySize,10);

differenceRates = cell(memorySize,1);
plusMinusRates = cell(memorySize,1);
timesNan = NaN(30000,12);
deltaNan = NaN(30000,12);

% Iterate through all subject directories in the chosen directory
fprintf("Beginning to process files for the %s task.\n",taskName);
for i=startDir:length(subjectDirectories)
    fprintf('Directory No. %i: %s\n',i, subjectDirectories(i).name);
    
    zipToLoad = dir2(fullfile(subjectDirectories(i).folder,subjectDirectories(i).name,fileToLoad));
    
    % Iterate through all trial files in the subject directory
    for j = 1:length(zipToLoad)
        data = zip_load(fullfile(zipToLoad(j).folder,zipToLoad(j).name));   % Loads the named file into a new structure called 'data'.
        data = KINARM_add_hand_kinematics(data);                            % Add hand velocity, acceleration and commanded forces to the data structure
        
        % Uncomment if you want to inspect which files are being loaded.
        % Usually this is too much output.
        % fprintf('File to load (%s): ''%s''\n',zipToLoad(j).name,data.file_label);
        
        % Asssess the trial file for various errors and edge cases
        if isempty(data)
            fprintf("This file produced an empty data struct! Skipping it and moving on to the next file.\n");
            movefile(fullfile(zipToLoad(j).folder,zipToLoad(j).name),fullfile(troubleRecordsFolder,zipToLoad(j).name));
            continue
        end
        
        if ~contains(data.file_label,["Object Hit","Object hit","object hit"])
            if verbose
                fprintf("This file is not an Object Hit task! Moving it to trouble records and moving on to the next file.\n");
            end
            movefile(fullfile(zipToLoad(j).folder,zipToLoad(j).name),fullfile(troubleRecordsFolder,zipToLoad(j).name));
            continue
        elseif OHflag && contains(data.file_label,["avoid","Avoid"])
            if verbose
                
                fprintf("This file is an Object-Hit-and-Avoid task! We were looking for Object-Hit tasks only.\n");
            end
            %movefile(fullfile(zipToLoad(j).folder,zipToLoad(j).name),fullfile(troubleRecordsFolder,zipToLoad(j).name));
            continue
        elseif ~OHflag && ~contains(data.file_label,["avoid","Avoid"])
            if verbose
                fprintf("This file isn't an Object-Hit-and-Avoid task! We are only looking for Object-Hit-and-Avoid tasks only right now.\n\n");
            end
            continue
        elseif turboFlag && ~contains(data.file_label,["turbo","Turbo","fast","Fast"])
            if verbose
                fprintf("This file isn't a Turbo variant of a task! We are only looking for Turbo variants right now.\n\n");
            end
            continue
        elseif ~turboFlag && contains(data.file_label,["turbo","Turbo"])
            if verbose
                fprintf("This file is a Turbo variant of a task! We are only looking for Classic variants right now.\n\n");
            end
            continue
        end
        
        % Process the trial file
        [decisionList,motorMetaData,steadyStateRate,overwhelmedHitRateMedian,overwhelmedHitRateIQR,steadyStateFrame,atLimitFrame,...
            decisionFrames,avgGap,targetCreates,~,hits,misses,drops,crashes,subjectAccuracy,earlyBinAccuracy,overwhelmedBinAccuracy,...
            earlyBinAvoidance,overwhelmedBinAvoidance] = processSingleOHFile(data,graphFlag,fileSaveType,verbose);
        
        % As long as we get real input back, save it!
        if ~isempty(decisionList)
            filesSeen = filesSeen + 1;
            dirFilesSeen = dirFilesSeen + 1;
            
            avgTargetCreationRate = movmean(targetCreates,1000,'omitnan')*200;
            avgHitRate = movmean(hits,1000,'omitnan')*200;
            avgSuccessRate = movmean(hits + drops,1000,'omitnan')*200;
            avgFailureRate = movmean(misses + crashes,1000,'omitnan')*200;
            
            % Calculate differential rate and median values
            differenceRate = avgTargetCreationRate - avgHitRate;
            plusMinusRate = avgSuccessRate - avgFailureRate;
            
            subjectNames(filesSeen) = subjectDirectories(i).name;
            trialNames(filesSeen) = regexp(zipToLoad(j).name,'\d\d\d\d-\d\d-\d\d_\d\d-\d\d-\d\d','match');
            differenceRates{filesSeen} = differenceRate;
            plusMinusRates{filesSeen} = plusMinusRate;
            newTimes = (1:length(differenceRate))/200;
            timesNan(:,filesSeen) = padarray(newTimes',30000-length(newTimes),NaN,'post');
            deltaNan(:,filesSeen) = padarray(differenceRate,30000-length(differenceRate),NaN,'post');
            decisionLists(filesSeen) = mat2cell(decisionList,size(decisionList,1),size(decisionList,2));
            decisionFramesList(filesSeen) = mat2cell(decisionFrames,size(decisionFrames,1),size(decisionFrames,2));
            metaDataArray(filesSeen) = copyMetaData(metaDataArray(filesSeen),motorMetaData);
            hitsScores(filesSeen) = sum(hits);
            dropsScores(filesSeen) = sum(drops);
            steadyStateRates(filesSeen) = steadyStateRate;
            owhrMedian(filesSeen) = overwhelmedHitRateMedian;
            owhrIQR(filesSeen) = overwhelmedHitRateIQR;
            steadyStateFrames(filesSeen) = steadyStateFrame;
            atLimitFrames(filesSeen) = atLimitFrame;
            avgGaps(filesSeen) = avgGap;
            maxAvgTargetCreationRate(filesSeen) = max(avgTargetCreationRate);
            subjectAccuracies(filesSeen) = subjectAccuracy;
            earlyBinAccuracies(filesSeen,:) = earlyBinAccuracy;
            overwhelmedBinAccuracies(filesSeen,:) = overwhelmedBinAccuracy;
            earlyBinAvoidances(filesSeen,:) = earlyBinAvoidance;
            overwhelmedBinAvoidances(filesSeen,:) = overwhelmedBinAvoidance;
        % If we don't get any real input back, note it and move on. 
        else
            fprintf("\nDecision List was empty! The task had to be skipped.\n");
        end
    end
    
    if dirFilesSeen > 0
        totalSubjects = totalSubjects + 1;
        dirFilesSeen = 0;
        fprintf("Processed %i trials from %i subjects.\n",filesSeen,totalSubjects);
    end
    
    if filesSeen > maxNumRecords
        blockFlag = true;
        break
    end
end

endDir = i;

if i < length(subjectDirectories)
    endFlag = 0;
else
    endFlag = 1;
end

trialSubjectString = sprintf("Directories %i to %i complete. Analysed a total of %i trials across %i different subjects.\n",startDir,endDir,filesSeen,totalSubjects);
fprintf(trialSubjectString);

% Trim arrays
subjectNames = subjectNames(1:filesSeen);
trialNames = trialNames(1:filesSeen);
decisionLists = decisionLists(1:filesSeen);
decisionFramesList = decisionFramesList(1:filesSeen);
metaDataArray = metaDataArray(1:filesSeen);
hitsScores = hitsScores(1:filesSeen);
dropsScores = dropsScores(1:filesSeen);
steadyStateRates = steadyStateRates(1:filesSeen);
owhrMedian = owhrMedian(1:filesSeen);
owhrIQR = owhrIQR(1:filesSeen);
steadyStateFrames = steadyStateFrames(1:filesSeen);
atLimitFrames = atLimitFrames(1:filesSeen);
avgGaps = avgGaps(1:filesSeen);
maxAvgTargetCreationRate = maxAvgTargetCreationRate(1:filesSeen);
subjectAccuracies = subjectAccuracies(1:filesSeen);
earlyBinAccuracies = earlyBinAccuracies(1:filesSeen,:);
overwhelmedBinAccuracies = overwhelmedBinAccuracies(1:filesSeen,:);
earlyBinAvoidances = earlyBinAvoidances(1:filesSeen,:);
overwhelmedBinAvoidances = overwhelmedBinAvoidances(1:filesSeen,:);

if blockFlag
    newEndString = sprintf('(%i-%i).mat',startDir,endDir);
    matFileName = strrep(matFileName,'.mat',newEndString);
end

% Save all of the results from the current block
save(matFileName,'subjectNames','trialNames','decisionLists','metaDataArray','hitsScores','dropsScores','steadyStateRates',...
    'owhrMedian','owhrIQR','steadyStateFrames','atLimitFrames','decisionFramesList','avgGaps','trialSubjectString',...
    'OHflag','turboFlag','maxAvgTargetCreationRate','subjectAccuracies','earlyBinAccuracies','overwhelmedBinAccuracies',...
    'earlyBinAvoidances','overwhelmedBinAvoidances');
fprintf('Finished saving data to the file ''%s''.\n',matFileName);

% If we wanted a diary, produce it
if diaryFlag
    fprintf('Beginning diary...\n');
    diary diaryName %#ok<*UNRCH>
    for i=1:length(decisionLists{1})
        d = decisionLists{1}{i};
        fprintf('%i. Decision: %i hits %i. Paddles: #1 [%i,%i] #2 [%i,%i].',decisionFramesList{1}(i),d{1},d{2},d{3},d{4});
        for j = 5:length(d)
            b = d{j};
            fprintf(' Ball %i [%i,%i] @ %i (%i).',b.id,b.xBin,b.yBin,b.speed,b.value);
        end
        fprintf('\n');
    end
    diary
    fprintf('...ending diary.\n');
end

% Only produce this global summary if we aren't analysing a block
if (graphFlag && ~blockFlag)
    fprintf('Beginning graph...\n');
    differenceRates = differenceRates(1:filesSeen);
    
    f = figure;
    colspace = linspace(0.5,0.9,filesSeen);
    meanTimes = nanmean(timesNan,2);
    xTimes = meanTimes(~isnan(meanTimes));
    meanDeltas = nanmean(deltaNan,2);
    yDeltas = meanDeltas(~isnan(meanDeltas));
    p = polyfit(xTimes,yDeltas,4);
    deltaFit = polyval(p,xTimes);
    
    hold on;
    titleString = sprintf('The difference between target creation and contact rates\nstays near zero until it doesn''t for %s',taskName);
    title(titleString);
    for i = 1:filesSeen
        rates = differenceRates{i};
        times = (1:length(rates))/200;
        if i==1 || i==filesSeen
            plot(times,rates,'Color',[0 colspace(i) colspace(i)]);
        else
            plot(times,rates,'Color',[0 colspace(i) colspace(i)],'HandleVisibility','off');
        end
    end
    plot(xTimes,deltaFit,'r','LineWidth', 2);
    xlabel("Time (s)");
    ylabel({'Difference in target creation','and contact rates (Hz)'});
    yline(0,'k','LineWidth',2,'HandleVisibility','off');
    nString = sprintf('(n = %i)',length(differenceRates));
    legend('Individual Trials',nString,'Average Difference','location','northwest');
    ax = gca;
    ax.FontSize = 24;
    f.WindowState = 'maximized';
    if filesFlag == 3 || filesFlag == 5 || filesFlag == 6 || filesFlag == 7
        if blockFlag
            fileSaveName = sprintf('BCCRatesNearZero(%s)(%i-%i).png',taskName,startDir,endDir);
        else
            fileSaveName = sprintf('BCCRatesNearZero(%s).png',taskName);
        end
        saveas(f,fileSaveName);
    end
    hold off;
    fprintf('...finished graph.\n');
end
end

%% FUNCTIONS

% Set a metadata structure for object-hit
function metaDataArray = setObjectHitMetaData()
metaDataArray.DISTRACTOR_HITS_LEFT = -1;
metaDataArray.DISTRACTOR_HITS_RIGHT = -1;
metaDataArray.DISTRACTOR_HITS_TOTAL = -1;
metaDataArray.HAND_BIAS_OF_HITS = -1;
metaDataArray.HAND_SELECTION_OVERLAP = -1;
metaDataArray.HAND_TRANSITION = -1;
metaDataArray.HAS_DISTRACTORS = 'false';
metaDataArray.HITS_LEFT = -1;
metaDataArray.HITS_RIGHT = -1;
metaDataArray.HITS_TOTAL = -1;
metaDataArray.MEAN_HAND_SPEED_BIAS = -1;
metaDataArray.MEAN_HAND_SPEED_LEFT = -1;
metaDataArray.MEAN_HAND_SPEED_RIGHT = -1;
metaDataArray.MEDIAN_ERROR = -1;
metaDataArray.MISS_BIAS = -1;
metaDataArray.MOVEMENT_AREA_BIAS = -1;
metaDataArray.MOVEMENT_AREA_LEFT = -1;
metaDataArray.MOVEMENT_AREA_RIGHT = -1;
metaDataArray.SCREEN_BOUNDS = [0;0;0;0];
metaDataArray.TOTAL_DISTRACTORS = -1;
metaDataArray.TOTAL_TARGETS = -1;
metaDataArray.DESCRIPTIONS = cell(1,21);
end

% Set a metadata structure for object-hit-and-avoid
function metaDataArray = setObjectHitAndAvoidMetaData()
metaDataArray.DISTRACTOR_HITS_LEFT = -1;
metaDataArray.DISTRACTOR_HITS_RIGHT = -1;
metaDataArray.DISTRACTOR_HITS_TOTAL = -1;
metaDataArray.DISTRACTOR_PROPORTION = -1; %%%
metaDataArray.HAND_BIAS_OF_HITS = -1;
metaDataArray.HAND_SELECTION_OVERLAP = -1;
metaDataArray.HAND_TRANSITION = -1;
metaDataArray.HAS_DISTRACTORS = 'false';
metaDataArray.HITS_LEFT = -1;
metaDataArray.HITS_RIGHT = -1;
metaDataArray.HITS_TOTAL = -1;
metaDataArray.MEAN_HAND_SPEED_BIAS = -1;
metaDataArray.MEAN_HAND_SPEED_LEFT = -1;
metaDataArray.MEAN_HAND_SPEED_RIGHT = -1;
metaDataArray.MEDIAN_ERROR = -1;
metaDataArray.MISS_BIAS = -1;
metaDataArray.MOVEMENT_AREA_BIAS = -1;
metaDataArray.MOVEMENT_AREA_LEFT = -1;
metaDataArray.MOVEMENT_AREA_RIGHT = -1;
metaDataArray.OBJECTS_HIT = -1; %%%
metaDataArray.OBJECT_PROCESSING_RATE = -1; %%%
metaDataArray.SCREEN_BOUNDS = [0;0;0;0];
metaDataArray.TOTAL_DISTRACTORS = -1;
metaDataArray.TOTAL_TARGETS = -1;
metaDataArray.DESCRIPTIONS = cell(1,21);
end
