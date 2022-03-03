%% This script produces subfigures/panels for the paper "Capacity Limits Lead to Information Bottlenecks for Ongoing Rapid Motor Behaviour"
%  Also produced are subfigures/panels for the associated technical report,
%  "Assessing Human Performance in Ongoing Rapid Motor Behaviours Tasks"
%
%  In some cases, additional plots are created for all tasks or for all
%  comparisons where only one is used in the paper.
%
%  2 December, 2021
%
%  Requires: Statistics and Machine Learning Toolbox
%
%#ok<*UNRCH>
%#ok<*DEFNU>

% File level settings
fileSaveType = 'png';
addpath('../Utilities/');
addpath('../Data Analysis/');

% Figure-specific flags allow a subset of subfigures/panels to be produced
figure4Flag = 0;
figure57Flag = 0;
figure56Flag = 0;
techRptFlag = 1;

%% Figure 4: Predictions made by simulations
if figure4Flag
    % Get access to the SSS Experiment Framework
    addpath('../SSS Simulation/');
    system_types = {'AW','MW','FC'};
    system_type_names = {'Additive Workload','Multiplicative Workload','Flexible Capacity'};
    
    % Load and do some book-keeping for data
    load('sss_results.mat');
    numTasks = size(tasks,1);           % The number of tasks included in sss_results
    baselineTaskNum = numTasks;         % The baseline task is always stored at the end
    minDifficulty = min(tasks,[],2);    % Calculate each task's minimum difficulty
    maxDifficulty = max(tasks,[],2);    % Calculate each task's maximum difficulty
    [sortDiff,minDifficultySort] = sortrows([minDifficulty maxDifficulty],[1 2]); % Sort tasks from most to least difficult
    
    capacityA = capacities(:,:,1);      % The capacities for all Module A's.
    capacityB = capacities(:,:,2);      % The capacities for all Module B's.
    
    % Demonstrate that all the models achieve a Steady State Rate
    % Additive Workload Model
    [~,input_record,output_record] = AdditiveWorkloadSerialSystem([3 2.5],[1 1],0);
    f = getSteadyStateFigure(input_record,output_record);
    set(groot,'currentfigure',f);
    hold on;
    yline(2.5,'k--','LineWidth',8);
    yline(3,'k--','LineWidth',8);
    hold off;
    savefigas(f,'SteadyStateRateSimulationBaselineTaskAW',fileSaveType);
    
    [~,input_record,output_record] = AdditiveWorkloadSerialSystem([3 2.5],[0.4 1],0);
    f = getSteadyStateFigure(input_record,output_record);
    set(groot,'currentfigure',f);
    hold on;
    yline(2.5,'k--','LineWidth',8);
    yline(1.2,'k--','LineWidth',8);
    hold off;
    savefigas(f,'SteadyStateRateSimulationSecondTaskAW',fileSaveType);
    
    % Multiplicative Workload Model
    [~,input_record,output_record] = MultiplicativeWorkloadSerialSystem([3 2.5],[1 1],0);
    f = getSteadyStateFigure(input_record,output_record);
    set(groot,'currentfigure',f);
    hold on;
    yline(2.5,'k--','LineWidth',8);
    yline(3,'k--','LineWidth',8);
    hold off;
    savefigas(f,'SteadyStateRateSimulationBaselineTaskMW',fileSaveType);
    
    [~,input_record,output_record] = MultiplicativeWorkloadSerialSystem([3 2.5],[0.4 1],0);
    f = getSteadyStateFigure(input_record,output_record);
    set(groot,'currentfigure',f);
    hold on;
    yline(2.5,'k--','LineWidth',8);
    yline(1.2,'k--','LineWidth',8);
    hold off;
    savefigas(f,'SteadyStateRateSimulationSecondTaskMW',fileSaveType);
    
    % Flexible Capacity Model
    [~,input_record,output_record] = FlexibleCapacitySerialSystem([3 2.5],[1 1],0);
    f = getSteadyStateFigure(input_record,output_record);
    set(groot,'currentfigure',f);
    hold on;
    yline(2.5,'k--','LineWidth',8);
    yline(3,'k--','LineWidth',8);
    hold off;
    savefigas(f,'SteadyStateRateSimulationBaselineTaskFC',fileSaveType);
    
    [~,input_record,output_record] = FlexibleCapacitySerialSystem([3 2.5],[0.4 1],0);
    f = getSteadyStateFigure(input_record,output_record);
    set(groot,'currentfigure',f);
    hold on;
    yline(2.5,'k--','LineWidth',8);
    yline(1.2,'k--','LineWidth',8);
    hold off;
    savefigas(f,'SteadyStateRateSimulationSecondTaskFC',fileSaveType);
    
    fprintf('We are considering the following tasks: Task 1 (%0.2f %0.2f); Task 2 (%0.2f %0.2f); and Task 3 (%0.2f %0.2f).\n',...
        tasks(baselineTaskNum,:),tasks(40,:),tasks(17,:));
    
    % Also assess and plot population-level mean performance of systems 
    % across all tasks for the Technical Report
    f = figure;
    hold on;
    
    % Count how many systems there are, since not all systems were
    % simulated.
    numSystems = sum(~isnan(results(:,:,1,1)),'all');
    
    % Calculate the mean performance and the standard error of the mean for
    % all Multiplicative Workload systems for each task.
    mwMeans = squeeze(mean(results(:,:,:,1),[1 2],'omitnan'));
    mwSEM = squeeze(std(results(:,:,:,1),0,[1 2],'omitnan'))./sqrt(numSystems);
    
    % Calculate the mean performance and the standard error of the mean for
    % all Flexible Capacity systems for each task.
    fcMeans = squeeze(mean(results(:,:,:,3),[1 2],'omitnan'));
    fcSEM = squeeze(std(results(:,:,:,3),0,[1 2],'omitnan'))./sqrt(numSystems);
    
    % Normalize performance and plot error bars.
    maxPerf = fcMeans(100);
    colour1 = [0.1059 0.6196 0.4667];
    colour2 = [0.4588 0.4392 0.7020];
    e1 = errorbar(mwMeans(minDifficultySort)/maxPerf,mwSEM,'LineStyle','none','LineWidth',10,'CapSize',15,'Color',[colour1 0.6]);
    e3 = errorbar(fcMeans(minDifficultySort)/maxPerf,fcSEM,'LineStyle','none','LineWidth',10,'CapSize',15,'Color',[colour2 0.6]);
    set([e1.Bar, e1.Line], 'ColorType', 'truecoloralpha', 'ColorData', [e1.Line.ColorData(1:3); 153]);
    set([e3.Bar, e3.Line], 'ColorType', 'truecoloralpha', 'ColorData', [e3.Line.ColorData(1:3); 153]);
   
    % Uncomment if you want to hightlight the tasks used for simple
    % comparisons.
    %xline(17,'--k','LineWidth',3,'HandleVisibility','off');
    %xline(40,'--k','LineWidth',3,'HandleVisibility','off');
    %xline(100,'--k','LineWidth',3,'HandleVisibility','off');

    ax = gca;
    ax.FontSize = 72;
    ax.FontWeight = 'bold';
    ax.LineWidth = 15;
    f.WindowState = 'maximized';
    hold off;
    savefigas(f,'AllSystemsMeanPerformanceAllTasks',fileSaveType);
    
    % For each model, produce a sample regression between the baseline task
    % and a second task. Then produce a summary of the regression
    % parameters for all comparisons between the baseline task and another
    % task.
    for i = 1:length(system_types)
        fprintf('Now analysing the %s type of parallel-serial system.\n',system_type_names{i});        
        fprintf('The mean performace level on Task 1 is %0.2f Hz.\n',mean(results(:,:,baselineTaskNum,i),'all','omitnan'));
        fprintf('The mean performace level on Task 2 is %0.2f Hz.\n',mean(results(:,:,40,i),'all','omitnan'));
        
        % Scatter plot comparison between Baseline Task and Task 2
        task1Results = results(:,:,baselineTaskNum,i);
        task2Results = results(:,:,40,i);
        [f,mdl1,mdl2] = scatterWithRegression(task1Results(:),task2Results(:),[0 3.5],[0 3.5],1);
        fprintf('The adjusted R^{2} value for the Task 1 v Task 2 comparison is %0.2f. The slope is %0.2f and the intercept is %0.2f.\n',mdl1.Rsquared.Adjusted,mdl1.Coefficients.Estimate(2),mdl1.Coefficients.Estimate(1));
        set(groot,'currentfigure',f);
        hold on;
        axis square;
        yticks([0 1 2 3])
        xticks([0 1 2 3])
        hold off;
        savefigas(f,sprintf('SystemOutputBetweenTasks%s',system_types{i}),fileSaveType);
        
        % Repeat this regression comparison between the Baseline Task and
        % all other Tasks, recording the Slope, Y-Intercept, and Adjusted
        % R-squared values for each regression.
        modelCoefficients = NaN(size(tasks,1),3);
        x = squeeze(results(:,:,baselineTaskNum,i));
        for j=1:size(tasks,1)
            y = squeeze(results(:,:,j,i));
            mdl = fitlm(x(:),y(:));
            modelCoefficients(j,1:2) = mdl.Coefficients.Estimate;
            modelCoefficients(j,3) = mdl.Rsquared.Adjusted;
        end
        
        % Plot the Regression Slope, Y-Intercept, and Adjusted R-squared
        % values for each task to demonstrate the effect of task difficulty
        % on these parameters.
        f = figure;
        hold on;
        colour1 = [0.1059 0.6196 0.4667];
        colour2 = [0.8510 0.3725 0.0078];
        colour3 = [0.4588 0.4392 0.7020];
        scatter(1:size(tasks,1),modelCoefficients(minDifficultySort,1),125,colour1,'filled');
        scatter(1:size(tasks,1),modelCoefficients(minDifficultySort,2),125,colour2,'filled');
        scatter(1:size(tasks,1),modelCoefficients(minDifficultySort,3),125,colour3,'filled');
        yline(1,'k','LineWidth',1,'HandleVisibility','off');
        yline(0,'k','LineWidth',6,'HandleVisibility','off');
        yline(-1,'k','LineWidth',1,'HandleVisibility','off');
        
        % Uncomment if you want to highlight the comparison tasks
        %xline(40,'--k','LineWidth',3,'HandleVisibility','off');
        %xline(100,'--k','LineWidth',3,'HandleVisibility','off');
        
        % Include the regions of practical equivalence for the Slope and
        % Y-intercept regression terms
        xCoords = [0 0 100 100];
        yCoords = [-0.12 0.12 0.12 -0.12];
        yCoords2 = [0.96 1.06 1.06 0.96];
        patch(xCoords,yCoords,[0.1 0.1 0.1],'FaceAlpha',0.3,'LineStyle','none','HandleVisibility','off');
        patch(xCoords,yCoords2,[0.1 0.1 0.1],'FaceAlpha',0.3,'LineStyle','none','HandleVisibility','off');
        
        ylim([-1.2 1.2]);
        yticks([-1 0 1]);
        ax = gca;
        ax.FontSize = 72;
        ax.FontWeight = 'bold';
        ax.LineWidth = 15;
        f.WindowState = 'maximized';
        hold off;
        savefigas(f,sprintf('taskComparisonCoefficients%s',system_types{i}),fileSaveType);
        
        % Confirm that everything matches up with the first scatter plot
        % comparison between the Baseline Task and Task 2.
        fprintf('The adjusted R^{2} value for the [1.0 1.0] v [0.7 1.0] comparison is %0.2f. The slope is %0.2f and the intercept is %0.2f.\n',modelCoefficients(40,3),modelCoefficients(40,2),modelCoefficients(40,1));        
    end
end

%% Figures 5-7: Comparing performance between OH and OHA-Comparing performance between tasks
if figure57Flag
    % In order to match naming conventions from other files in this
    % project, we append '-SingleTrial' to the file name to denote that the
    % data depicted contains only the first pair of trials from each
    % participant.
    singleTrialString = '-SingleTrial';
    
    % Prepare references for all four comparisons; only the first three are
    % included in the paper.
    taskOnes = {'OHCombinedValues.mat','decisionListsCompOH.mat','turboOHCombinedValues.mat','decisionListsCompTOH.mat'};
    taskTwos = {'OHACombinedValues.mat','decisionListsCompTOH.mat','turboOHACombinedValues.mat','OHACombinedValues.mat'};
    comparisons = {'OH-OHA','OH-TOH','TOH-TOHA','TOH-OHA'};
    
    % For each comparison, 
    for z=1:length(taskOnes)
        
        % Load the data for the current comparison.
        t1Results = load(taskOnes{z});
        t2Results = load(taskTwos{z});
        
        % Set long and short task names based on the loaded data.
        if t1Results.OHflag
            t1Name = 'Object-Hit';
            t1ShortName = 'OH';
        else
            t1Name = 'Object-Hit-and-Avoid';
            t1ShortName = 'OHA';
        end
        if t2Results.OHflag
            t2Name = 'Object-Hit';
            t2ShortName = 'OH';
        else
            t2Name = 'Object-Hit-and-Avoid';
            t2ShortName = 'OHA';
        end
        if t1Results.turboFlag
            t1Name = ['Turbo ' t1Name];
            t1ShortName = ['T' t1ShortName];
        end
        if t2Results.turboFlag
            t2Name = ['Turbo ' t2Name];
            t2ShortName = ['T' t2ShortName];
        end
        
        fprintf('\nComparison is %s vs. %s.\n', t1ShortName,t2ShortName);
        
        % Determine which subjects in trial 1 have also performed trial 2.
        % From there, determine based on trial time-stamps how to align the
        % trials.
        t1UniqueSubjects = unique(t1Results.subjectNames);  % The unique subject names from Task 1
        estLength = length(t1Results.steadyStateRates);     % An estimated length for creating data structures
        trial1SSRs = zeros(estLength,1);                    % Record the SSRs for the trials from Task 1
        trial2SSRs = zeros(estLength,1);                    % Record the SSRs for the trials from Task 2
        trial1Max = zeros(estLength,1);                     % Record the maximum average target creation rate for the trials from Task 1
        trial2Max = zeros(estLength,1);                     % Record the maximum average target creation rate for the trials from Task 2
        matchedSubjects = strings(estLength,1);             % Record the subjects who were successfully matched between Tasks
        idx = 1;                                            % A counter pointing to the next empty slot in all of the above data structures
        
        % For each nuique subject
        for i=1:length(t1UniqueSubjects)
            % Get their name
            subjectName = t1UniqueSubjects(i);
            
            % Find their trials
            t1TrialIndices = find(ismember(t1Results.subjectNames,subjectName));
            t2TrialIndices = find(ismember(t2Results.subjectNames,subjectName));
            
            % Determine if the day of the two trials should match for
            % matching trials. 
            if z == 2 || z == 4
                dayMatchFlag = 0;
            else
                dayMatchFlag = 1;
            end
            
            % Get all of the algined trial indices for Task 1 and Task 2
            [t1TrialIndices,t2TrialIndices] = alignTrials(t1TrialIndices,t2TrialIndices,t1Results.trialNames,t2Results.trialNames,0,dayMatchFlag,0);
            
            % If both returned arrays are non-empty, then record the
            % results
            if ~isempty(t1TrialIndices) && ~isempty(t2TrialIndices)
                trial1SSRs(idx) = t1Results.steadyStateRates(t1TrialIndices(1));
                trial2SSRs(idx) = t2Results.steadyStateRates(t2TrialIndices(1));
                trial1Max(idx) = t1Results.maxAvgTargetCreationRate(t1TrialIndices(1));
                trial2Max(idx) = t2Results.maxAvgTargetCreationRate(t2TrialIndices(1));
                matchedSubjects(idx) = subjectName;
                idx = idx + 1;
            end
        end
        
        % Trim the arrays as necessary
        trial1SSRs = trial1SSRs(1:idx-1);
        trial2SSRs = trial2SSRs(1:idx-1);
        trial1Max = trial1Max(1:idx-1);
        trial2Max = trial2Max(1:idx-1);
        matchedSubjects = matchedSubjects(1:idx-1);
        
        % Get the peak target creation rate across both trial types
        maxTargetCreationRate = max([trial1Max;trial2Max]);
        
        % Get the minimum SSRs across both trial types
        minSSR = min([trial1SSRs;trial2SSRs]);
        maxSSR = max([trial1SSRs;trial2SSRs]);
        
        % Set axis limits for the figure
        xLimits = [0 4];
        if (max(trial2SSRs) < 2)
            yLimits = [0 2];
        else
            yLimits = [0 3];
        end
        
        if strcmp(t1ShortName,'OH')
            t1ROPE = 0.30;
        elseif strcmp(t1ShortName,'TOH')
            t1ROPE = 0.48;
        else
            fprintf('Unrecognized Task 1 Name! Expected either OH or TOH but got %s instead.\n',t1ShortName);
        end
        
        [f,mdl] = scatterWithRegression(trial1SSRs,trial2SSRs,xLimits,yLimits,1);
        set(groot,'currentfigure',f);
        hold on;
        minCoord = min([xLimits(1)  yLimits(1)]);
        maxCoord = max([xLimits(2) yLimits(2)]);
        
        % Draw the region of practical equivalence around the unity line
        ropeDist = sqrt((t1ROPE^2)/2);
        xCoords = [(minCoord+ropeDist) (minCoord-ropeDist) (maxCoord-ropeDist) (maxCoord+ropeDist)];
        yCoords = [(minCoord-ropeDist) (minCoord+ropeDist) (maxCoord+ropeDist) (maxCoord-ropeDist)];
        patch(xCoords,yCoords,[0.1 0.1 0.1],'FaceAlpha',0.3,'LineStyle','none');
        
        xlim(xLimits);
        ylim(yLimits);
        hold off;
        savefigas(f,sprintf('CompareSSRs(%s-%s)',t1ShortName,t2ShortName),fileSaveType);
        
        fprintf('The adjusted R2 value for SSRs between %s is %0.2f. The intercept is %0.2f and the slope is %0.2f.\n',...
            comparisons{z},mdl.Rsquared.Adjusted,mdl.Coefficients.Estimate(1),mdl.Coefficients.Estimate(2));
        fprintf('This comparison had %i trials in it.\n',length(trial1SSRs));
        
        % Save steady-state rates as a csv file for Bayesian Data Analysis in R
        fileNameString = sprintf('SubjectSteadyStateRates(%s-%s%s).csv',t1ShortName,t2ShortName,singleTrialString);
        csvwrite(fileNameString,[trial1SSRs trial2SSRs]);
    end
end

%% Figure 5-6: Drops in performance between tasks are especially significant for better performers
if figure56Flag
    % Set-up some values that have been calculated elsewhere in the project
    ohROPE = -0.3;
    tohROPE = -0.48;
    ohSSRsd = 0.425001;
    tohSSRsd = 0.369793;
    
    % Create probability distribution functions for the OH ROPE, OH 95% CI
    % for healthy subjects, TOH ROPE, and TOH 95% CI for healthy subjects.
    ohROPEpdf = makedist('Uniform',ohROPE,-ohROPE);
    oh95HealthyPDF = makedist('Normal',0,ohSSRsd);
    tohROPEpdf = makedist('Uniform',tohROPE,-tohROPE);
    toh95HealthyPDF = makedist('Normal',0,tohSSRsd);
    
    % Plot population-level PDF comparisons
    % Mean SSR differences for each comparison
    ohohaPopMean = -0.7334;
    tohtohaPopMean = -0.5316;
    ohtohPopMean = -0.5305;
    % Create probability distribution functions for each comparison
    ohohaDiff = makedist('Normal','mu',-0.7334,'sigma',0.2434);
    tohtohaDiff = makedist('Normal','mu',-0.5316,'sigma',0.2571);
    ohtohDiff = makedist('Normal','mu',-0.5305,'sigma',0.2544);
    
    % Create and save population-level PDF comparisons
    f = plotPDFsWithPopulationMean(ohROPEpdf,oh95HealthyPDF,ohohaDiff);
    savefigas(f,'ohohaPopPDFs',fileSaveType);
    f = plotPDFsWithPopulationMean(tohROPEpdf,toh95HealthyPDF,tohtohaDiff);
    savefigas(f,'tohtohaPopPDFs',fileSaveType);
    f = plotPDFsWithPopulationMean(ohROPEpdf,oh95HealthyPDF,ohtohDiff);
    savefigas(f,'ohtohPopPDFs',fileSaveType);
    
    fprintf('The OH-OHA posterior distribution is %0.2f%% below the ROPE. The posterior mean difference has a z-score of %0.4f.\n',...
        100*cdf(ohohaDiff,ohROPE),ohohaPopMean/ohSSRsd);
    fprintf('The TOH-TOHA posterior distribution is %0.2f%% below the ROPE. The posterior mean difference has a z-score of %0.4f.\n',...
        100*cdf(tohtohaDiff,tohROPE),tohtohaPopMean/tohSSRsd);
    fprintf('The OH-TOH posterior distribution is %0.2f%% below the ROPE. The posterior mean difference has a z-score of %0.4f.\n',...
        100*cdf(ohtohDiff,ohROPE),ohtohPopMean/ohSSRsd);
    
    % Plot quartile-level comparisons
    % OH-OHA comparison
    ohohaQuartileMus = [-0.5078,-0.6825,-0.7828,-0.9607];
    ohohaQuartileSigmas = [0.1786 0.1608 0.1913 0.1931];
    ohohaDiffQ1 = makedist('Normal','mu',ohohaQuartileMus(1),'sigma',ohohaQuartileSigmas(1));
    ohohaDiffQ2 = makedist('Normal','mu',ohohaQuartileMus(2),'sigma',ohohaQuartileSigmas(2));
    ohohaDiffQ3 = makedist('Normal','mu',ohohaQuartileMus(3),'sigma',ohohaQuartileSigmas(3));
    ohohaDiffQ4 = makedist('Normal','mu',ohohaQuartileMus(4),'sigma',ohohaQuartileSigmas(4));
    
    f = plotPDFsWithPopulationMean(ohROPEpdf,oh95HealthyPDF,ohohaDiffQ1,ohohaPopMean);
    savefigas(f,'ohohaQ1PDFs',fileSaveType);
    f = plotPDFsWithPopulationMean(ohROPEpdf,oh95HealthyPDF,ohohaDiffQ2,ohohaPopMean);
    savefigas(f,'ohohaQ2PDFs',fileSaveType);
    f = plotPDFsWithPopulationMean(ohROPEpdf,oh95HealthyPDF,ohohaDiffQ3,ohohaPopMean);
    savefigas(f,'ohohaQ3PDFs',fileSaveType);
    f = plotPDFsWithPopulationMean(ohROPEpdf,oh95HealthyPDF,ohohaDiffQ4,ohohaPopMean);
    savefigas(f,'ohohaQ4PDFs',fileSaveType);
    
    fprintf('The OH-OHA Quartile 1 posterior distribution is %0.2f%% below the ROPE. The posterior mean difference has a z-score of %0.4f.\n',...
        100*cdf(ohohaDiffQ1,ohROPE),ohohaQuartileMus(1)/ohSSRsd);
    fprintf('The OH-OHA Quartile 2 posterior distribution is %0.2f%% below the ROPE. The posterior mean difference has a z-score of %0.4f.\n',...
        100*cdf(ohohaDiffQ2,ohROPE),ohohaQuartileMus(2)/ohSSRsd);
    fprintf('The OH-OHA Quartile 3 posterior distribution is %0.2f%% below the ROPE. The posterior mean difference has a z-score of %0.4f.\n',...
        100*cdf(ohohaDiffQ3,ohROPE),ohohaQuartileMus(3)/ohSSRsd);
    fprintf('The OH-OHA Quartile 4 posterior distribution is %0.2f%% below the ROPE. The posterior mean difference has a z-score of %0.4f.\n',...
        100*cdf(ohohaDiffQ4,ohROPE),ohohaQuartileMus(4)/ohSSRsd);
    
    % TOH-TOHA comparison
    tohtohaQuartileMeans = [-0.3585 -0.4500 -0.5412 -0.7738];
    tohtohaQuartileSigmas = [0.1638 0.1846 0.1989 0.2641];
    tohtohaDiffQ1 = makedist('Normal','mu',tohtohaQuartileMeans(1),'sigma',tohtohaQuartileSigmas(1));
    tohtohaDiffQ2 = makedist('Normal','mu',tohtohaQuartileMeans(2),'sigma',tohtohaQuartileSigmas(2));
    tohtohaDiffQ3 = makedist('Normal','mu',tohtohaQuartileMeans(3),'sigma',tohtohaQuartileSigmas(3));
    tohtohaDiffQ4 = makedist('Normal','mu',tohtohaQuartileMeans(4),'sigma',tohtohaQuartileSigmas(4));
    
    f = plotPDFsWithPopulationMean(tohROPEpdf,toh95HealthyPDF,tohtohaDiffQ1,tohtohaPopMean);
    savefigas(f,'tohtohaQ1PDFs',fileSaveType);
    f = plotPDFsWithPopulationMean(tohROPEpdf,toh95HealthyPDF,tohtohaDiffQ2,tohtohaPopMean);
    savefigas(f,'tohtohaQ2PDFs',fileSaveType);
    f = plotPDFsWithPopulationMean(tohROPEpdf,toh95HealthyPDF,tohtohaDiffQ3,tohtohaPopMean);
    savefigas(f,'tohtohaQ3PDFs',fileSaveType);
    f = plotPDFsWithPopulationMean(tohROPEpdf,toh95HealthyPDF,tohtohaDiffQ4,tohtohaPopMean);
    savefigas(f,'tohtohaQ4PDFs',fileSaveType);
    
    fprintf('The TOH-TOHA Quartile 1 posterior distribution is %0.2f%% below the ROPE. The posterior mean difference has a z-score of %0.4f.\n',...
        100*cdf(tohtohaDiffQ1,tohROPE),tohtohaQuartileMeans(1)/tohSSRsd);
    fprintf('The TOH-TOHA Quartile 2 posterior distribution is %0.2f%% below the ROPE. The posterior mean difference has a z-score of %0.4f.\n',...
        100*cdf(tohtohaDiffQ2,tohROPE),tohtohaQuartileMeans(2)/tohSSRsd);
    fprintf('The TOH-TOHA Quartile 3 posterior distribution is %0.2f%% below the ROPE. The posterior mean difference has a z-score of %0.4f.\n',...
        100*cdf(tohtohaDiffQ3,tohROPE),tohtohaQuartileMeans(3)/tohSSRsd);
    fprintf('The TOH-TOHA Quartile 4 posterior distribution is %0.2f%% below the ROPE. The posterior mean difference has a z-score of %0.4f.\n',...
        100*cdf(tohtohaDiffQ4,tohROPE),tohtohaQuartileMeans(4)/tohSSRsd);
    
    % OH-TOH comparison
    ohtohQuartileMeans = [-0.4211,-0.4270,-0.5179,-0.7442];
    ohtohQuartileSigmas = [0.1898 0.2076 0.1908 0.3690];
    ohtohDiffQ1 = makedist('Normal','mu',ohtohQuartileMeans(1),'sigma',ohtohQuartileSigmas(1));
    ohtohDiffQ2 = makedist('Normal','mu',ohtohQuartileMeans(2),'sigma',ohtohQuartileSigmas(2));
    ohtohDiffQ3 = makedist('Normal','mu',ohtohQuartileMeans(3),'sigma',ohtohQuartileSigmas(3));
    ohtohDiffQ4 = makedist('Normal','mu',ohtohQuartileMeans(4),'sigma',ohtohQuartileSigmas(4));
    
    f = plotPDFsWithPopulationMean(ohROPEpdf,oh95HealthyPDF,ohtohDiffQ1,ohtohPopMean);
    savefigas(f,'ohtohQ1PDFs',fileSaveType);
    f = plotPDFsWithPopulationMean(ohROPEpdf,oh95HealthyPDF,ohtohDiffQ2,ohtohPopMean);
    savefigas(f,'ohtohQ2PDFs',fileSaveType);
    f = plotPDFsWithPopulationMean(ohROPEpdf,oh95HealthyPDF,ohtohDiffQ3,ohtohPopMean);
    savefigas(f,'ohtohQ3PDFs',fileSaveType);
    f = plotPDFsWithPopulationMean(ohROPEpdf,oh95HealthyPDF,ohtohDiffQ4,ohtohPopMean);
    savefigas(f,'ohtohQ4PDFs',fileSaveType);
    
    fprintf('The OH-TOH Quartile 1 posterior distribution is %0.2f%% below the ROPE. The posterior mean difference has a z-score of %0.4f.\n',...
        100*cdf(ohtohDiffQ1,ohROPE),ohtohQuartileMeans(1)/ohSSRsd);
    fprintf('The OH-TOH Quartile 2 posterior distribution is %0.2f%% below the ROPE. The posterior mean difference has a z-score of %0.4f.\n',...
        100*cdf(ohtohDiffQ2,ohROPE),ohtohQuartileMeans(2)/ohSSRsd);
    fprintf('The OH-TOH Quartile 3 posterior distribution is %0.2f%% below the ROPE. The posterior mean difference has a z-score of %0.4f.\n',...
        100*cdf(ohtohDiffQ3,ohROPE),ohtohQuartileMeans(3)/ohSSRsd);
    fprintf('The OH-TOH Quartile 4 posterior distribution is %0.2f%% below the ROPE. The posterior mean difference has a z-score of %0.4f.\n',...
        100*cdf(ohtohDiffQ4,ohROPE),ohtohQuartileMeans(4)/ohSSRsd);
    
end

%% Technical Report
if techRptFlag
    % To maintain naming consistency with other files in this project, we
    % will append '-AllTrials' to the end of figure names to indicate that
    % the figures are based on all available trials for a given task.
    singleTrialString = '-AllTrials';
    
    % We do not want to produce our scatter plots with the  unity line pre-built
    unityFlag = 0;
    
    % Object-Hit
    taskName = 'Object-Hit';
    taskShortName = 'OH';
    
    temp1 = load('OHCombinedValues.mat','hitsScores','dropsScores','steadyStateRates','metaDataArray','earlyBinAccuracies','overwhelmedBinAccuracies');
    temp2 = load('decisionListsCompOH.mat','hitsScores','dropsScores','steadyStateRates','metaDataArray','earlyBinAccuracies','overwhelmedBinAccuracies');
    data = mergeStructs(temp1,temp2);
    
    % The range of targets hit by subjects across all trials
    f = rangeOfHitsHistogram(data.hitsScores,taskName);
    fileNameString = sprintf('RangeOfHits(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    fprintf('Number of target hits for OH: Range: %i - %i, Median %i, Mean %0.2f. Total trials: %i.\n',min(data.hitsScores),max(data.hitsScores),median(data.hitsScores),mean(data.hitsScores),length(data.hitsScores));

    % The per-bin accuracy of subjects during their early and overwhelmed phases
    f = binAccuracyPlots(data.earlyBinAccuracies,data.overwhelmedBinAccuracies);
    fileNameString = sprintf('BinAccuracies(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    
    % Comparing steady-state rates with the number of targets hit
    f = scatterWithRegression(data.hitsScores,data.steadyStateRates,[100 300],[0 4],unityFlag);
    fileNameString = sprintf('HitsSSRScatter(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    
    % Object-Hit-and-Avoid
    taskName = 'Object-Hit-and-Avoid';
    taskShortName = 'OHA';
    
    data = load('OHACombinedValues.mat','hitsScores','dropsScores','steadyStateRates',...
        'metaDataArray','earlyBinAccuracies','overwhelmedBinAccuracies','earlyBinAvoidances','overwhelmedBinAvoidances');
    
    % The range of targets hit by subjects across all trials 
    f = rangeOfHitsHistogram(data.hitsScores,taskName);
    fileNameString = sprintf('RangeOfHits(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    fprintf('Number of target hits for OHA: Range: %i - %i, Median %i, Mean %0.2f. Total trials: %i.\n',min(data.hitsScores),max(data.hitsScores),median(data.hitsScores),mean(data.hitsScores),length(data.hitsScores));
    

    % The per-bin accuracy of subjects during their early and overwhelmed phases
    f = binAccuracyPlots(data.earlyBinAccuracies,data.overwhelmedBinAccuracies);
    fileNameString = sprintf('BinAccuracies(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    
    % Comparing steady-state rates with the number of targets hit
    f = scatterWithRegression(data.hitsScores,data.steadyStateRates,[0 200],[0 4],unityFlag);
    fileNameString = sprintf('HitsSSRScatter(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    
    % Comparing subjects' success at hitting targets and avoiding
    % distractors
    xLimits = [50 200];
    yLimits = [30 100];
    [f,mdl] = scatterWithRegression(data.hitsScores,data.dropsScores,xLimits,yLimits,2);
    fprintf('The adjusted R2 value for hits vs. drops in OHA is %0.2f. The intercept is %0.2f and the slope is %0.2f.\n',mdl.Rsquared.Adjusted,mdl.Coefficients.Estimate(1),mdl.Coefficients.Estimate(2));
    fileNameString = sprintf('TargetDistractorCorrelation(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    
    
    % Turbo Object-Hit
    taskName = 'Turbo Object-Hit';
    taskShortName = 'TOH';
    
    temp1 = load('turboOHCombinedValues.mat','hitsScores','dropsScores','steadyStateRates','metaDataArray','earlyBinAccuracies','overwhelmedBinAccuracies');
    temp2 = load('decisionListsCompTOH.mat','hitsScores','dropsScores','steadyStateRates','metaDataArray','earlyBinAccuracies','overwhelmedBinAccuracies');
    data = mergeStructs(temp1,temp2);
    
    % The range of targets hit by subjects across all trials
    f = rangeOfHitsHistogram(data.hitsScores,taskName);
    fileNameString = sprintf('RangeOfHits(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    fprintf('Number of target hits for TOH: Range: %i - %i, Median %i, Mean %0.2f. Total trials: %i.\n',min(data.hitsScores),max(data.hitsScores),median(data.hitsScores),mean(data.hitsScores),length(data.hitsScores));
    
    % The per-bin accuracy of subjects during their early and overwhelmed phases
    f = binAccuracyPlots(data.earlyBinAccuracies,data.overwhelmedBinAccuracies);
    fileNameString = sprintf('BinAccuracies(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    
    % Comparing steady-state rates with the number of targets hit
    f = scatterWithRegression(data.hitsScores,data.steadyStateRates,[100 300],[0 4],unityFlag);
    fileNameString = sprintf('HitsSSRScatter(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    
    % Turbo Object-Hit-and-Avoid
    taskName = 'Turbo Object-Hit-and-Avoid';
    taskShortName = 'TOHA';
    xLabel = sprintf('Number of Targets Hit (n)');
    yLabel = sprintf('Number of Distractors Avoided (n)');
    
    data = load('turboOHACombinedValues.mat','hitsScores','dropsScores','steadyStateRates','metaDataArray','earlyBinAccuracies','overwhelmedBinAccuracies');
    
    % The range of targets hit by subjects across all trials
    f = rangeOfHitsHistogram(data.hitsScores,taskName);
    fileNameString = sprintf('RangeOfHits(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    fprintf('Number of target hits for TOHA: Range: %i - %i, Median %i, Mean %0.2f. Total trials: %i.\n',min(data.hitsScores),max(data.hitsScores),median(data.hitsScores),mean(data.hitsScores),length(data.hitsScores));
    
    % The per-bin accuracy of subjects during their early and overwhelmed phases
    f = binAccuracyPlots(data.earlyBinAccuracies,data.overwhelmedBinAccuracies);
    fileNameString = sprintf('BinAccuracies(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    
    % Comparing steady-state rates with the number of targets hit
    f = scatterWithRegression(data.hitsScores,data.steadyStateRates,[0 200],[0 4],unityFlag);
    fileNameString = sprintf('HitsSSRScatter(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    
    % Comparing subjects' success at hitting targets and avoiding
    % distractors
    xLimits = [50 200];
    yLimits = [30 100];
    [f,mdl] = scatterWithRegression(data.hitsScores,data.dropsScores,xLimits,yLimits,2);
    fprintf('The adjusted R2 value for hits vs. drops in TOHA is %0.2f. The intercept is %0.2f and the slope is %0.2f.\n',mdl.Rsquared.Adjusted,mdl.Coefficients.Estimate(1),mdl.Coefficients.Estimate(2));
    fileNameString = sprintf('TargetDistractorCorrelation(%s%s)',taskShortName,singleTrialString);
    savefigas(f,fileNameString,fileSaveType);
    
    % Produce an example Steady-State Rate Figure
    % Another file in this project is capable of producing these for all of
    % the trial we have data for, but that seems excessive for right now.
    exampleFile = '/home/richard/Documents/School/PhD/Scott Lab/KINARM data/Classic OH-OHA (control)/OH/ACB_ACB_3725/3725_2017-03-31_11-34-21.zip';
    numObjects = 300;
    
    data = zip_load(exampleFile);                       % Loads the named file into a new structure called 'data'
    data = KINARM_add_hand_kinematics(data);            % Add hand velocity, acceleration and commanded forces to the data structure
    [eventsList,~] = getEventsList(data,numObjects);    % Extract the germane events from the trial
    [targetCreates,distractorCreates,hits,misses,drops,crashes] = getBasicStatistics(data,eventsList); % Computes basic statistics from the trial
    
    % We multiply the event-occurence arrays by 200 to transform their averages
    % from 'per-frame' to 'per-second'.
    targetCreates = targetCreates * 200;
    hits = hits * 200;
    
    f = getSteadyStateFigure(targetCreates,hits);
    savefigas(f,'SteadyStateRateExample',fileSaveType);
end

% rangeOfHitsHistogram is used to produce histograms of the number of
% targets hit by subjects for various tasks. It knows that different types
% of tasks have different numbers of targets available and trials recorded,
% so it adapts accordingly.
function f = rangeOfHitsHistogram(hitsScores,taskName)
switch taskName
    case 'Object-Hit'
        maxHits = 300;
        maxSubj = 150;
    case 'Object-Hit-and-Avoid'
        maxHits = 200;
        maxSubj = 150;
    case 'Turbo Object-Hit'
        maxHits = 300;
        maxSubj = 700;
    case 'Turbo Object-Hit-and-Avoid'
        maxHits = 200;
        maxSubj = 700;
    otherwise
        fprintf('Unexpected task name! Got ''%s''.\n',taskName);
        keyboard
end

hitsRange = 0:(maxHits/20):maxHits;

f = figure;
hold on;
histogram(hitsScores,hitsRange,'FaceColor',[0 0 0],'FaceAlpha',0.6,'EdgeAlpha',0.6)
axis([0 maxHits 0 maxSubj]);
ax = gca;
ax.FontSize = 72;
ax.FontWeight = 'bold';
ax.LineWidth = 10;
axis square
f.WindowState = 'maximized';
hold off;
end

% binAccuracyPlots plots the mean per-bin accuracy of all subjects in their
% respective early and overwhelmed phases. It also computes and plots the
% standard errors of these means, but in practice these are not visible.
function f = binAccuracyPlots(earlyBinAccuracies,overwhelmedBinAccuracies)

binAccuracies = NaN(10,2);
binSEM = NaN(10,2);
for i=1:size(binAccuracies,1)
    binAccuracies(i,1) = mean(earlyBinAccuracies(:,i),'omitNan');
    binSEM(i,1) = std(earlyBinAccuracies(:,i),'omitnan')/sqrt(sum(~isnan(earlyBinAccuracies(:,i))));
    binAccuracies(i,2) = mean(overwhelmedBinAccuracies(:,i),'omitNan');
    binSEM(i,2) = std(overwhelmedBinAccuracies(:,i),'omitnan')/sqrt(sum(~isnan(overwhelmedBinAccuracies(:,i))));
end

f = figure;
hold on;
colour1 = [0.0742 0.3867 0.6172];
colour2 = [0.6172 0.1445 0.0742];
errorbar(1:10,binAccuracies(:,1),binSEM(:,1),'LineStyle','none','LineWidth',20,'CapSize',30,'Marker','o','MarkerSize',20,'MarkerFaceColor','auto','Color',[colour1 0.6]);
errorbar(1:10,binAccuracies(:,2),binSEM(:,2),'LineStyle','none','LineWidth',20,'CapSize',30,'Marker','o','MarkerSize',20,'MarkerFaceColor','auto','Color',[colour2 0.6]);
xlim([0.5 10.5]);
ylim([0 1]);
xticks([1 4 7 10]);
yticks([0 1]);
ax = gca;
ax.FontSize = 72;
ax.FontWeight = 'bold';
ax.LineWidth = 10;
f.WindowState = 'maximized';
hold off;
end

% scatterWithRegression produces a straightforward scatter plot and also
% includes the linear regression line. If unityFlag is set, then it also
% draws the unity line for comparison.
function [f,mdl] = scatterWithRegression(x,y,xLimits,yLimits,unityFlag)
f = figure;
hold on;

% Scatter
scatter(x,y,150,'black','filled','MarkerFaceAlpha',0.6);
pause(1);

% Unity Line
switch unityFlag
    case 1
        minCoord = min([xLimits(1)  yLimits(1)]);
        maxCoord = max([xLimits(2) yLimits(2)]);
        plot([minCoord maxCoord],[minCoord maxCoord],'k','LineWidth',6)
    case 2
        plot([0 200],[0 100],'k','LineWidth',6);
end

% Regression
mdl = fitlm(x,y,'linear');
xVals = linspace(xLimits(1),xLimits(2));
[mdlVals,~] = predict(mdl,xVals','Alpha',0.05);

plot(xVals,mdlVals,'Color',[0.8500 0.3250 0.0980],'LineWidth',10);

xlim(xLimits);
ylim(yLimits);

ax = gca;
ax.FontSize = 72;
ax.FontWeight = 'bold';
ax.LineWidth = 15;
f.WindowState = 'maximized';
hold off;
end

% getSteadyStateFigure smooths and then plots the target creation and
% target hit rates for an individual trial.
function f = getSteadyStateFigure(targetCreates,hits)
% Implement Kolmogorov-Zurbenko filters: 2 iterations, window sizes 1001
% data points
avgTargetCreationRate = movmean(movmean(targetCreates,1001,'omitnan'),1001,'omitnan');
avgHitRate = movmean(movmean(hits,1001,'omitnan'),1001,'omitnan');

% Calculate where the average object creation rate begins to fall off, this
% is the "end of the task" where the subject is no longer being
% increasingly stressed.
[~,endOfTask] = max(avgTargetCreationRate);

% Calculate differential rate and median values
differenceRate = avgTargetCreationRate(1:endOfTask) - avgHitRate(1:endOfTask);
medianDifferenceRate = movmedian(differenceRate,1001);
madDifferenceRate = movmad(differenceRate,1001);

% Determine the last frame for which the lower bound was below 0.1
steadyStateFrame = 1;
for i=endOfTask:-1:1
    if (medianDifferenceRate(i) - madDifferenceRate(i)) < 0.1
        steadyStateFrame = i;
        break
    end
end

% Calculate the steady state rate
steadyTimes = (steadyStateFrame:endOfTask)/200;
overwhelmedSteadyStateRate = mean(hits(steadyStateFrame:endOfTask));

% Catch some edge cases. If the steady state rate is NaN, that's a
% problem. If the steady state rate is 0 then the subject never got
% overwhelmed and we use the maximum hit rate instead.
if isnan(overwhelmedSteadyStateRate)
    minkeyboard
elseif (length(steadyTimes) < 200) || (overwhelmedSteadyStateRate == 0)
    overwhelmedSteadyStateRate = max(avgHitRate);
    steadyStateFrame = endOfTask;
end

% Keep track of all candidate steady state frames for graphing results
candidateSteadyStateFrames = steadyStateFrame;

% The Overwhelmed Phase: the steadyStateRate is the number of balls contacted
while 1
    % Check to see if overwhelmed phase should start earlier based on the
    % subject performing at a higher level than the calculated steady state
    % rate calculated.
    candidateFrame = steadyStateFrame;
    for i=steadyStateFrame:-1:1
        if (avgHitRate(i) > mean(hits(i:endOfTask)))
            candidateFrame = i;
            candidateSteadyStateRate = mean(hits(i:endOfTask));
            
            % If we are within 5 seconds of starting, give the search a chance
            % chance to keep finding "better" windows.
        elseif i < (steadyStateFrame - 1000)
            break
        end
        
    end
    
    % Record the current candidate frame if it is different
    if candidateSteadyStateFrames(end) ~= candidateFrame
        candidateSteadyStateFrames = [candidateSteadyStateFrames candidateFrame]; %#ok<*AGROW>
    end
    
    % If the subject was overwhelmed at this new candidate frame, keep
    % pulling back to see where they weren't overhelmed.
    newStartFrame = candidateFrame;
    for i=newStartFrame:-1:1
        if (medianDifferenceRate(i) - madDifferenceRate(i)) > 0.1
            candidateFrame = i;
            candidateSteadyStateRate = mean(hits(candidateFrame:endOfTask));
            
            % If we are within 5 seconds of starting, give the search a chance
            % chance to keep finding "better" windows.
        elseif i < (newStartFrame - 1000)
            break
        end
    end
    
    % Record the current candidate frame if it is different
    if candidateSteadyStateFrames(end) ~= candidateFrame
        candidateSteadyStateFrames = [candidateSteadyStateFrames candidateFrame];
    end
    
    % Update steady state frame if using the candidate frame rewards the
    % subject with a higher steady state rate.
    if candidateFrame < steadyStateFrame
        steadyStateFrame = candidateFrame;
        steadyTimes = (steadyStateFrame:endOfTask)/200;
        overwhelmedSteadyStateRate = candidateSteadyStateRate;
    else
        % Calculate comparison steady state values
        steadyStateRate = overwhelmedSteadyStateRate;
        break
    end
end

x = (1:endOfTask)/200;
f = figure('WindowState','maximized');
set(gcf,'Visible','on');
hold on;
plot(x,avgTargetCreationRate(1:endOfTask),'Color','black','LineWidth',10,'DisplayName','Target Creation Rate');
plot(x,avgHitRate(1:endOfTask),'Color',[106 142 210]/255,'LineWidth',10,'DisplayName','Target Contact Rate');
plot(steadyTimes,steadyStateRate*ones(length(steadyTimes),1),'Color',[165 24 24]/255,'LineWidth',8,'DisplayName','Overwhelmed Phase');

ax = gca;
ax.FontSize = 72;
ax.FontWeight = 'bold';
ax.LineWidth = 15;
hold off;
end

% plotPDFsWithPopulationMean produces figures with (potentially
% overlapping) probability distribution functions with their respective
% means highlighted.
function f = plotPDFsWithPopulationMean(rope,healthy95,posterior,varargin)
x = -1.5:0.05:1.5;

f = figure;
hold on;
area(x,pdf(rope,x),'FaceColor',"#58d68d",'FaceAlpha',0.6,'EdgeColor',"#58d68d");
area(x,pdf(healthy95,x),'FaceColor',"#495ab0",'FaceAlpha',0.6,'EdgeColor',"#495ab0");
area(x,pdf(posterior,x),'FaceColor',"#e84133",'FaceAlpha',0.6,'EdgeColor',"#e84133");
xline(mean(posterior),'Color',"#e84133",'LineWidth',10);
if nargin == 4
    xline(varargin{1},'LineWidth',10);
end

xlim([-1.5 1.5]);
ylim([0 2.5]);
yticks([]);
xticks([-1 0 1]);

ax = gca;
ax.FontSize = 48;
ax.FontWeight = 'bold';
ax.LineWidth = 10;
f.WindowState = 'maximized';
hold off;
end