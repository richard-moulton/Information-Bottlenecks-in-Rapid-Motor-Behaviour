%% A script for combining the desired values from multiple .mat files.
%
%  3 December, 2021

function matFileCombiner(varargin)
if nargin == 4
    ohListNames = varargin{1};
    ohaListNames = varargin{2};
    tohListNames = varargin{3};
    tohaListNames = varargin{4};
% Allows running the function as a script for the homebrew set-up.
else
    ohListNames = {'fullDecisionListsOH(1-384).mat','fullDecisionListsOH(385-620).mat'};
    ohaListNames = {'fullDecisionListsOHA(1-376).mat','fullDecisionListsOHA(377-516).mat'};
    tohListNames = {'fullDecisionListsTurboOH(1-245).mat','fullDecisionListsTurboOH(246-490).mat','fullDecisionListsTurboOH(491-747).mat',...
        'fullDecisionListsTurboOH(748-990).mat','fullDecisionListsTurboOH(991-1233).mat','fullDecisionListsTurboOH(1234-1498).mat',...
        'fullDecisionListsTurboOH(1499-1669).mat'};
    tohaListNames = {'fullDecisionListsTurboOHA(1-246).mat','fullDecisionListsTurboOHA(247-494).mat','fullDecisionListsTurboOHA(495-747).mat',...
        'fullDecisionListsTurboOHA(748-990).mat','fullDecisionListsTurboOHA(991-1233).mat','fullDecisionListsTurboOHA(1234-1499).mat',...
        'fullDecisionListsTurboOHA(1500-1669).mat'};
end

allListNames = {ohListNames,ohaListNames,tohListNames,tohaListNames};
saveFileNames = {'OHCombinedValues.mat','OHACombinedValues.mat','turboOHCombinedValues.mat','turboOHACombinedValues.mat'};

for i=1:length(allListNames)
    fprintf("\nTask %i: Processing for %s...\n",i,saveFileNames{i});
    
    % Initialize arrays
    subjectNames = strings(1,1);
    trialNames = strings(1,1);
    hitsScores = [];
    dropsScores = [];
    steadyStateRates = [];
    maxAvgTargetCreationRate = [];
    metaDataArray = [];
    earlyBinAccuracies = [];
    overwhelmedBinAccuracies = [];
    earlyBinAvoidances = [];
    overwhelmedBinAvoidances = [];
    
    % For every .mat file, append its data to the new arrays.
    if ~isempty(allListNames{i})
        for j=1:length(allListNames{i})
            fprintf("%i: Batch %s\n",j,allListNames{i}{j});
            temp = load(allListNames{i}{j},'subjectNames','trialNames','hitsScores',...
                'dropsScores','steadyStateRates','maxAvgTargetCreationRate',...
                'OHflag','turboFlag','metaDataArray','earlyBinAccuracies',...
                'overwhelmedBinAccuracies','earlyBinAvoidances','overwhelmedBinAvoidances');
            
            subjectNames = [subjectNames; temp.subjectNames];
            trialNames = [trialNames; temp.trialNames];
            hitsScores = [hitsScores; temp.hitsScores];
            dropsScores = [dropsScores; temp.dropsScores];
            steadyStateRates = [steadyStateRates; temp.steadyStateRates]; %#ok<*AGROW>
            maxAvgTargetCreationRate = [maxAvgTargetCreationRate; temp.maxAvgTargetCreationRate];
            OHflag = temp.OHflag;
            turboFlag = temp.turboFlag;
            metaDataArray = [metaDataArray,temp.metaDataArray];
            earlyBinAccuracies = [earlyBinAccuracies; temp.earlyBinAccuracies];
            overwhelmedBinAccuracies = [overwhelmedBinAccuracies; temp.overwhelmedBinAccuracies];
            earlyBinAvoidances = [earlyBinAvoidances; temp.earlyBinAvoidances];
            overwhelmedBinAvoidances = [overwhelmedBinAvoidances; temp.overwhelmedBinAvoidances];
        end
        
        % Get rid of the blank first entry for strings
        subjectNames = subjectNames(2:end);
        trialNames = trialNames(2:end);
        
        % Save the new .mat file.
        save(saveFileNames{i},'subjectNames','trialNames','hitsScores','dropsScores',...
            'steadyStateRates','maxAvgTargetCreationRate',...
            'OHflag','turboFlag','metaDataArray','earlyBinAccuracies','overwhelmedBinAccuracies',...
            'earlyBinAvoidances','overwhelmedBinAvoidances');
        fprintf("Saved %s\n",saveFileNames{i});
    else
        fprintf('Did not save %s - the respective list of mat files was empty!\n',saveFileNames{i});
    end
end
end