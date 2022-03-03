%% This function aligns trials based on their date-time distance.
%  This is useful for ensuring that paired trials are being analysed (if
%  that is what is desired).
%
%  3 Decemeber 2021
%
%  %#ok<*ST2NM>

function [t1Matched, t2Matched] = alignTrials(t1Idx,t2Idx,t1TrialsFull,t2TrialsFull,singleTrialFlag,dayMatchFlag,verbose)

if verbose
    fprintf("\n");
end

% Pull out the specific trial names needed.
t1Trials = t1TrialsFull(t1Idx);
t2Trials = t2TrialsFull(t2Idx);
t1Matched = zeros(length(t1Idx),1);
t2Matched = zeros(length(t2Idx),1);
numMatched = 0;

if verbose
    fprintf("Trial 1:\n");
    for i=1:length(t1Trials)
        fprintf("%s\n",t1Trials(i));
    end
    
    fprintf("Trial 2:\n");
    for i=1:length(t2Trials)
        fprintf("%s\n",t2Trials(i));
    end
    
    fprintf("Matched Trials:\n");
end

% If dayMatchFlag is set, then we require two trials to occur on the same
% day before we match them. This is helpful to find those trials that occur
% back-to-back on the same day and drop out any "oddball" trials that are
% recorded.
if dayMatchFlag
    if length(t1Trials) <= length(t2Trials)
        for i=1:length(t1Trials)
            [date,time] = getDateTime(t1Trials(i),verbose);
            
            [matchFlag,matchIdx] = checkDateTime(t2Trials,date,time,verbose);
            
            if matchFlag
                numMatched = numMatched + 1;
                t1Matched(numMatched) = t1Idx(i);
                t2Matched(numMatched) = t2Idx(matchIdx);
                if verbose
                    fprintf("%s - %s\n",t1Trials(i),t2Trials(matchIdx));
                end
            end
        end
    else
        for i=1:length(t2Trials)
            [date,time] = getDateTime(t2Trials(i),verbose);
            
            [matchFlag,matchIdx] = checkDateTime(t1Trials,date,time,verbose);
            
            if matchFlag
                numMatched = numMatched + 1;
                t2Matched(numMatched) = t2Idx(i);
                t1Matched(numMatched) = t1Idx(matchIdx);
                if verbose
                    fprintf("%s - %s\n",t1Trials(matchIdx),t2Trials(i));
                end
            end
        end
    end
% If dayMatchFlag is not set, then we drop the requirement for trials to
% occur on the same day. This is necessary for data sets where subjects
% have previously completed Task 1 and are then brought back to complete
% Task 2.
else
    numMatched = min(length(t1Trials),length(t2Trials));
    t1Matched(1:numMatched) = t1Idx(1:numMatched);
    t2Matched(1:numMatched) = t2Idx(1:numMatched);
    if verbose
        for z=1:numMatched
            fprintf("%s - %s\n",t1Trials(z),t2Trials(z));
        end
    end
end

% Trim the arrays as necessary
t1Matched = t1Matched(1:numMatched);
t2Matched = t2Matched(1:numMatched);

% If singleTrialFlag is set, we want only the first pair of trials.
if singleTrialFlag
    [t1Matched,t2Matched] = getEarliestPair(t1Matched,t2Matched,t1TrialsFull,verbose);
end
end

% Based on a trial string, return the date and time portions.
function [date,time] = getDateTime(trialName,verbose)
dateExp = '(?<year>\d{4})-(?<month>\d{2})-(?<day>\d{2})_';
date = regexp(trialName,dateExp,'names');
date.year = str2num(date.year);
date.month = str2num(date.month);
date.day = str2num(date.day);

timeExp = '_(?<hour>\d{2})-(?<minute>\d{2})-(?<second>\d{2})';
time = regexp(trialName,timeExp,'names');
time.hour = str2num(time.hour);
time.minute = str2num(time.minute);
time.second = str2num(time.second);

if verbose
    fprintf('%s gives %2d-%2d-%2d for the date and %2d:%2d:%2d for the time.\n',trialName,date.year,date.month,date.day,time.hour,time.minute,time.second);
end
end

% Check if a particular trial's date and time match well with the trial
% dates and times in an argument list. Day must match exactly, if there are
% multiple options then closest in time is used.
function [matchFlag,matchIdx] = checkDateTime(trialList,date,time,verbose)

numTrials = length(trialList);
matchDist = zeros(numTrials,1);
matchIdx = zeros(numTrials,1);
numMatches = 0;

for i=1:numTrials
    [cDate,cTime] = getDateTime(trialList(i),verbose);
    
    if((date.year == cDate.year) && (date.month == cDate.month) && (date.day == cDate.day))
        numMatches = numMatches + 1;
        matchDist(numMatches) = getTimeDistance(time,cTime,1,verbose);
        matchIdx(numMatches) = i;
    end
end

if numMatches == 0
    matchFlag = 0;
    matchIdx = [];
elseif numMatches == 1
    matchFlag = 1;
    matchIdx = matchIdx(1);
else
    matchFlag = 1;
    [~,matchIdx] = min(matchDist);
end
end

% A function that computes the distance in days between two date/time
% groups
function dateTimeDist = getDateTimeDistance(date1,time1,date2,time2,verbose)
dateTimeDist = 365*(date1.year-date2.year) + 30*(date1.month-date2.month) + date1.day-date2.day + getTimeDistance(time1,time2,0,verbose)/86400;
end

% A function that computes the number of seconds between two time groups
function timeDist = getTimeDistance(time1,time2,absoluteDistance,verbose)
timeDist = 3600*(time1.hour - time2.hour) + 60*(time1.minute - time2.minute) + (time1.second - time2.second);
if absoluteDistance
    timeDist = abs(timeDist);
end
if verbose
    fprintf('Distance between %i:%i:%i and %i:%i:%i was calculated to be %i seconds.\n', time1.hour,time1.minute,time1.second,time2.hour,time2.minute,time2.second,timeDist);
end
end

% A function that reduces a lists of matched trials to their earliest
% respective trials
function [matched1,matched2] = getEarliestPair(matched1,matched2,t1TrialNames,verbose)
if length(matched1) > 1
    [earliestDate,earliestTime] = getDateTime(t1TrialNames(matched1(1)),verbose);
    earliestIdx = 1;
    
    for i=2:length(matched1)
        [newDate,newTime] = getDateTime(t1TrialNames(matched1(i)),verbose);
        
        if verbose
            fprintf("Comparing Index %i (%i-%i-%i @ %i:%i:%i) and Index %i (%i-%i-%i @ %i:%i:%i)\n",earliestIdx,...
                earliestDate.year,earliestDate.month,earliestDate.day,earliestTime.hour,earliestTime.minute,earliestTime.second,...
                i,newDate.year,newDate.month,newDate.day,newTime.hour,newTime.minute,newTime.second);
        end
        
        if getDateTimeDistance(earliestDate,earliestTime,newDate,newTime,verbose) > 0
            earliestDate = newDate;
            earliestTime = newTime;
            earliestIdx = i;
        end
    end
    
    matched1 = matched1(earliestIdx);
    matched2 = matched2(earliestIdx);
end
end