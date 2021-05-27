clc;
clear all;
close all;

problem.cost_func = @(x) minOne(x);
problem.n_var = 100;

params.max_it = 500;
params.n_pop = 500;
params.percent_child = 1;
params.Mu = 0.07;
params.beta = 1;

out = binaryGeneticAlgorithm(problem, params);

figure;
plot(out.best_cost, 'LineWidth', 2);
xlabel('Iterations');
ylabel('Best Cost');
grid on;
