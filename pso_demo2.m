clc;
clear all;
close all;

problem.cost_func = @(x) minMax(x);
problem.n_var = 5;
problem.var_min = 0.1;
problem.var_max = 1;

params.max_it = 1000;
params.n_pop = 40;
params.w = 1;
params.w_damp = 0.99;
params.c1 = 2;
params.c2 = 2;
params.show_iter_info = true;

out = particleSwarm(problem, params);

pop = out.pop;
global_best = out.global_best;
best_costs = out.best_costs;

figure;
semilogy(best_costs, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Costs');
