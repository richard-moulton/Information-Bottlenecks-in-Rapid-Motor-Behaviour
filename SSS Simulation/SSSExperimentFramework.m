%% A script set up to run experiments using the SimpleSerialSystem set-up.
%
%  Specify the task environments based on the difficulty modifiers for the
%  different modules.
%  Specify the module processing capacities for the various systems.
%  Run simulations for all acceptable capacity settings and record the
%  results.
%
%  3 December, 2021

saveFlag = 1;               % Flag indicating whether or not to save the results
totalSimulations = 0;       % Counter for number of simulations run

system_types = {'Additive Workload','Multiplicative Workload','Flexible Capacity'};
types_of_systems = length(system_types);       % Number of system setups being simulated

% Difficulty modifiers for different tasks, resulting in 100 different
% task environments
difficulties = 0.1:0.1:1;
[t1,t2] = meshgrid(difficulties,difficulties);
c = cat(2,t1,t2);
tasks = reshape(c,[],2);

% Module processing capacities, resulting in 625 possible different
% systems. Not all of these will be used, however.
module_capacities = 1:0.1:3.5;
num_module_capacities = length(module_capacities);

% Initialize arrays to store experiment parameters and results
results = NaN(num_module_capacities,num_module_capacities,size(tasks,1),types_of_systems);
capacities = NaN(num_module_capacities,num_module_capacities,2);
minCapacities = NaN(num_module_capacities,num_module_capacities);
maxCapacities = NaN(num_module_capacities,num_module_capacities);

% Experiment-level noise parameter. Doesn't really change the results of
% the analysis, but avoids such straight/smooth lines in all the figures.
noiseLevel = 0.05;

% Iterate through all of the possible combinations for module capacities
for i = 1:num_module_capacities
    for j = 1:num_module_capacities
        
        % Only run the simulations if one module isn't going to throttle
        % the other's processing capacity.
        if abs(module_capacities(i) - module_capacities(j)) > (0.33 * module_capacities(i))
            fprintf('Simulation #XX. Set-up %i-%i-XX. Capacities %0.2f and %0.2f were too far apart.\n',i,j,module_capacities(i),module_capacities(j));
            continue
        else
            
            totalSimulations = totalSimulations + 1;
            
            % Again, this doesn't really change the results of the
            % analysis. Jittering the module capacities does prevent
            % systems from being plotted on top of each other in scatter
            % plots, though.
            final_module_capacities = [(module_capacities(i)+noiseLevel*randn(1)) (module_capacities(j)+noiseLevel*randn(1))];
            
            % Record the module capacities that have been generated.
            capacities(i,j,:) = final_module_capacities;
            minCapacities(i,j) = min(final_module_capacities);
            maxCapacities(i,j) = max(final_module_capacities);
            
            % Use the generated module capacities to run simulations in
            % every single task environment.
            for l = 1:size(tasks,1)
                fprintf('Simulation #%i. Set-up %i-%i-%i. System Capacities: %0.2f (A) and %0.2f (B). Module Mean Limits: %0.2f (A) and %0.2f (B). Task Difficulties: %0.2f (A) and %0.2f (B).\n',...
                    totalSimulations,i,j,l,module_capacities(i),module_capacities(j),final_module_capacities,tasks(l,:));
                visualizeFlag = 0;
                results(i,j,l,1) = AdditiveWorkloadSerialSystem(squeeze(capacities(i,j,:))',tasks(l,:),visualizeFlag);
                results(i,j,l,2) = MultiplicativeWorkloadSerialSystem(squeeze(capacities(i,j,:))',tasks(l,:),visualizeFlag);
                results(i,j,l,3) = FlexibleCapacitySerialSystem(squeeze(capacities(i,j,:))',tasks(l,:),visualizeFlag);
            end
        end
    end
end

% Save the experiment parameters and results.
if saveFlag
    save('sss_results.mat','results','capacities','minCapacities','maxCapacities','tasks','system_types');
end