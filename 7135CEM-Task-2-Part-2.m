clc;
clear;
close all;

% Define CEC'2005 Benchmark Functions
sphereFunction = @(x) sum(x.^2);
rosenbrockFunction = @(x) sum(100 * (x(2:end) - x(1:end-1).^2).^2 + (1 - x(1:end-1)).^2);
rastriginFunction = @(x) sum(x.^2 - 10 * cos(2 * pi * x) + 10);

% Problem Definition
dim = 10;  
lb = -5.12; 
ub = 5.12;  
numRuns = 15; 

functions = {sphereFunction, rosenbrockFunction, rastriginFunction};
functionNames = {'Sphere', 'Rosenbrock', 'Rastrigin'};

% Storage for results
resultsFmincon = zeros(numRuns, length(functions));
resultsPSO = zeros(numRuns, length(functions));
resultsSA = zeros(numRuns, length(functions));

for fIdx = 1:length(functions)
    objFunction = functions{fIdx};
    
    figure;
    hold on;
    title(['Optimization Convergence - ', functionNames{fIdx}]);
    xlabel('Iteration');
    ylabel('Objective Value');
    grid on;

    % Run fmincon 15 times
    for i = 1:numRuns
        options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off');
        x0 = rand(1, dim) * (ub - lb) + lb;
        [~, fval, ~, output] = fmincon(objFunction, x0, [], [], [], [], lb * ones(1, dim), ub * ones(1, dim), [], options);
        resultsFmincon(i, fIdx) = fval;
        plot(1:output.iterations, output.iterations, 'b-', 'LineWidth', 1);
    end
    
    % Run Custom PSO 15 times
    for i = 1:numRuns
        [~, fval, iterHistory] = custom_pso(objFunction, dim, lb, ub);
        resultsPSO(i, fIdx) = fval;
        plot(1:length(iterHistory), iterHistory, 'r-', 'LineWidth', 1);
    end
    
    % Run Custom Simulated Annealing 15 times
    for i = 1:numRuns
        [~, fval, iterHistory] = custom_sa(objFunction, dim, lb, ub);
        resultsSA(i, fIdx) = fval;
        plot(1:length(iterHistory), iterHistory, 'g-', 'LineWidth', 1);
    end
    
    legend('Fmincon', 'PSO', 'SA');
    
    % Compute statistics
    statsFmincon = [mean(resultsFmincon(:, fIdx)), std(resultsFmincon(:, fIdx)), min(resultsFmincon(:, fIdx)), max(resultsFmincon(:, fIdx))];
    statsPSO = [mean(resultsPSO(:, fIdx)), std(resultsPSO(:, fIdx)), min(resultsPSO(:, fIdx)), max(resultsPSO(:, fIdx))];
    statsSA = [mean(resultsSA(:, fIdx)), std(resultsSA(:, fIdx)), min(resultsSA(:, fIdx)), max(resultsSA(:, fIdx))];
    
    % Display results
    fprintf('\nOptimization Results for %s Function:\n', functionNames{fIdx});
    fprintf('Algorithm\tMean\tStd Dev\tBest\tWorst\n');
    fprintf('Fmincon\t%.6f\t%.6f\t%.6f\t%.6f\n', statsFmincon);
    fprintf('PSO\t%.6f\t%.6f\t%.6f\t%.6f\n', statsPSO);
    fprintf('SA\t%.6f\t%.6f\t%.6f\t%.6f\n', statsSA);
end

% Box Plot for performance comparison
figure;
boxplot([resultsFmincon(:), resultsPSO(:), resultsSA(:)], {'Fmincon', 'PSO', 'SA'});
title('Comparison of Optimization Algorithms');
ylabel('Final Objective Value');
grid on;

% ---------------------------------
%        CUSTOM PSO FUNCTION
% ---------------------------------
function [bestPosition, bestFitness, iterHistory] = custom_pso(objFunc, dim, lb, ub)
    swarmSize = 30;
    maxIter = 100;
    w = 0.7;  
    c1 = 1.5;  
    c2 = 1.5;  
    iterHistory = zeros(1, maxIter);

    particles = rand(swarmSize, dim) * (ub - lb) + lb;
    velocities = zeros(swarmSize, dim);
    pBest = particles;
    pBestFitness = arrayfun(@(i) objFunc(particles(i, :)), 1:swarmSize);
    
    [bestFitness, bestIdx] = min(pBestFitness);
    gBest = pBest(bestIdx, :);

    for iter = 1:maxIter
        for i = 1:swarmSize
            velocities(i, :) = w * velocities(i, :) + ...
                              c1 * rand * (pBest(i, :) - particles(i, :)) + ...
                              c2 * rand * (gBest - particles(i, :));
                          
            particles(i, :) = particles(i, :) + velocities(i, :);
            particles(i, :) = max(min(particles(i, :), ub), lb);
            
            newFitness = objFunc(particles(i, :));
            if newFitness < pBestFitness(i)
                pBest(i, :) = particles(i, :);
                pBestFitness(i) = newFitness;
            end
        end
        
        [newBestFitness, newBestIdx] = min(pBestFitness);
        if newBestFitness < bestFitness
            bestFitness = newBestFitness;
            gBest = pBest(newBestIdx, :);
        end
        iterHistory(iter) = bestFitness;
    end

    bestPosition = gBest;
end

% ---------------------------------
%  CUSTOM SIMULATED ANNEALING (SA)
% ---------------------------------
function [bestPosition, bestFitness, iterHistory] = custom_sa(objFunc, dim, lb, ub)
    maxIter = 100;
    T = 1000;
    alpha = 0.95;
    iterHistory = zeros(1, maxIter);
    
    bestPosition = rand(1, dim) * (ub - lb) + lb;
    bestFitness = objFunc(bestPosition);
    currentPosition = bestPosition;
    currentFitness = bestFitness;

    for iter = 1:maxIter
        newPosition = currentPosition + (rand(1, dim) - 0.5) * (ub - lb) * 0.1;
        newPosition = max(min(newPosition, ub), lb);
        newFitness = objFunc(newPosition);
        
        if newFitness < currentFitness || rand < exp((currentFitness - newFitness) / T)
            currentPosition = newPosition;
            currentFitness = newFitness;
            if newFitness < bestFitness
                bestFitness = newFitness;
                bestPosition = newPosition;
            end
        end
        
        T = T * alpha;
        iterHistory(iter) = bestFitness;
    end
end