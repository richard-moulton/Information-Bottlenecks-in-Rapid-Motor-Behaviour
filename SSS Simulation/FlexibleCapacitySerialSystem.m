%% Flexible Capacity Serial System
%  A parallel-serial processing system where a part of the system is a
%  limited pool of central resources that can be flexibly assigned to
%  modules.
%
%  3 December, 2021

function [finalRate,input_record,output_record] = FlexibleCapacitySerialSystem(systemCapacities,taskModifiers,visualize)
simulationLength = 140; % in seconds
simulationFrameRate = 200; % in Hz (ie frames/second)

simulationFrameTotal = simulationLength*simulationFrameRate + 1;

% The input module, time varying input that starts at these values
input_function= @(t) 0.5 + (0.025 * t) + randn(1)/10;
input_record = NaN(simulationFrameTotal,1);

% The output module simply records all system outputs
output_record = NaN(simulationFrameTotal,1);

% Initialize
in = input_function(1/200);
input_record(1) = in;

% Reduce each module's processing capacity equally to produce a central
% resource pool. This ensures that all comparable systems have the same
% "total" processing ability.
systemCapacities = systemCapacities - (1/length(systemCapacities));
centralResourcePool = 1;

moduleCapacities = systemCapacities + 0.05*randn([1 length(systemCapacities)]);
taskDifficulties = taskModifiers + 0.01*randn([1 length(taskModifiers)]);

moduleFrameCapacity = getModuleFrameCapacities(moduleCapacities,taskDifficulties,centralResourcePool);

output_record(1) = min([in moduleFrameCapacity]);

stepVector = (1:99)/100;

% Run the simulation
for i = 1:(simulationFrameTotal/100)
    in = input_function(i/2);
    input_record(i*100 + 1) = in;
    input_record((i*100)-98:(i*100)) = input_record((i*100)-99) + (input_record((i*100)-99) - input_record((i*100)+1))*stepVector;
    
    moduleCapacities = systemCapacities + 0.05*randn([1 length(systemCapacities)]);
    taskDifficulties = taskModifiers + 0.05*randn([1 length(taskModifiers)]);
    moduleFrameCapacity = getModuleFrameCapacities(moduleCapacities,taskDifficulties,centralResourcePool);
    
    output_record(i*100 + 1) = max([0 min([(in + 0.2*in*randn(1)) moduleFrameCapacity])]);
    output_record((i*100)-98:(i*100)) = output_record((i*100)-99) + (output_record((i*100)-99) - output_record((i*100)+1))*stepVector;
    
    if visualize
        fprintf('Frame %i: Considering input %.2f and module limits %s. Output was %.2f.\n',(i*100)+1,in,num2str(moduleFrameCapacity),output_record((i*100)+1));
    end
end

time = (1:simulationFrameTotal)/simulationFrameRate;

% Implement Kolmogorov-Zurbenko filters: 2 iterations, window sizes 1001
% data points
smooth_input_record = movmean(movmean(input_record,1001,'omitnan'),1001,'omitnan');
smooth_output_record = movmean(movmean(output_record,1001,'omitnan'),1001,'omitnan');

finalRate = mean(smooth_output_record((simulationFrameTotal - 5000):simulationFrameTotal));

% If we want to visualize the results of the simulation, then we produce a
% "steady-state rate" type of figure for the system with the various module
% processing capacities illustrated as horizontal lines.
if visualize
    systemCapacities = systemCapacities + (1/length(systemCapacities));
    
    f = figure;
    hold on;
    plot(time,smooth_input_record,'-k','LineWidth',2,'DisplayName','System Input Rate');
    plot(time,smooth_output_record,'-b','LineWidth',2,'DisplayName','System Output Rate');
    
    zz = 1;
    moduleIDs = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    for mod = 1:length(systemCapacities)
        yline(systemCapacities(mod),'DisplayName',sprintf('Module %s Mean Processing Capacity (%0.2f)',moduleIDs(zz),systemCapacities(mod)));
        zz = zz + 1;
    end
    
    lgd = legend('location','Best');
    
    title(sprintf('Task (%s), Limited Capacity System (%s + 1)',num2str(taskModifiers),num2str(systemCapacities)));
    xlabel('Time (s)');
    ylabel('Input and Output (a.u.)');
    ylim([0 4]);
    
    ax = gca;
    ax.FontSize = 32;
    ax.FontWeight = 'bold';
    ax.LineWidth = 2.5;
    lgd.FontSize = 24;
    lgd.FontWeight = 'normal';
    f.WindowState = 'maximized';
    hold off;
end
end

% A function for sharing the central resource pool between the modules as 
% necessary throughout the simulation
function moduleFrameCapacity = getModuleFrameCapacities(moduleCapacities,taskDifficulties,centralResourcePool)
moduleFrameCapacity = moduleCapacities.*taskDifficulties;


while centralResourcePool > 0
    [sortedCapacities,sortedIndices] = sort(moduleFrameCapacity);
    
    firstMin = sortedCapacities(1);
    minIndex = sortedIndices(1);
    secondMin = sortedCapacities(2);
    
    if (secondMin - firstMin) > centralResourcePool
        moduleFrameCapacity(minIndex) = firstMin + centralResourcePool;
        centralResourcePool = 0;
    elseif secondMin == firstMin
        indexes = find(moduleFrameCapacity==firstMin);
        moduleFrameCapacity(indexes) = firstMin + centralResourcePool/length(indexes);
        centralResourcePool = 0;
    else
        moduleFrameCapacity(minIndex) = secondMin;
        centralResourcePool = centralResourcePool - (secondMin - firstMin);
    end
end
end