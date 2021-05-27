clc;
clear all;
close all;

problem.cost_func = @(x) minMax(x);
problem.n_var = 100;

params.max_it = 100;
params.n_pop = 100;
params.percent_child = 1;
params.Mu = 0.7;
params.beta = 1;

out = binaryGeneticAlgorithm(problem, params);

figure;
plot(out.best_cost, 'LineWidth', 2);
xlabel('Iterations');
ylabel('Best Cost');
grid on;
