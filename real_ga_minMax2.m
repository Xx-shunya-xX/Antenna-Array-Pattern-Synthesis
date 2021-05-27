clc;
clear all;
close all;

problem.cost_func = @(x) minMax(x);
problem.n_var = 5;
problem.var_min = 0.1;
problem.var_max = 1;

params.max_it = 100; %
params.n_pop = 40; %
params.percent_child = 1;
params.Mu = 0.02;
params.beta = 1;
params.sigma = 0.1;

out = realCodedGeneticAlgorithm(problem, params);

disp('Positions : ');
disp(out.global_best.position);
figure;
semilogy(out.best_cost, 'LineWidth', 2);
xlabel('Iterations');
ylabel('Best Cost');
grid on;
