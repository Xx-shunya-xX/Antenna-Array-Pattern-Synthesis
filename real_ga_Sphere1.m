clc;
clear all;
close all;

problem.cost_func = @(x) minMax(x);
problem.n_var = 5;
problem.var_min = 0;
problem.var_max = 1;

params.max_it = 1000; %
params.n_pop = 1000; %
params.percent_child = 1;
params.Mu = 0.02;
params.beta = 1;
params.sigma = 0.1;

out = realCodedGeneticAlgorithm(problem, params);

figure;
semilogy(out.best_cost, 'LineWidth', 2);
xlabel('Iterations');
ylabel('Best Cost');
grid on;
