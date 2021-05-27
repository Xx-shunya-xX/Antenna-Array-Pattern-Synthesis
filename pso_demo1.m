clc;
clear all;
close all;

%% Problem Definition
problem.cost_func = @(x) minMax(x);	% Cost Function
problem.n_var = 5;					% Number of Unknown (Decision) Variables
problem.var_min = 0.1;				% Lower Bound of Decision
problem.var_max = 1;				% Upper Bound of Decision

%% Parameters of PSO
params.max_it = 15;					% Maximum Number of Iterations
params.n_pop = 40;					% Population (Swarm) Size

w_max = 0.9;
w_min = 0.5;

params.w = w_min - (w_max - w_min) / params.max_it;						% Inertia Coefficient
params.w_damp = 0.99;				% Damping Ration of Inertia Weight
params.c1 = 2;						% Personal Acceleration Coefficient
params.c2 = 2;						% Social Acceleration Coefficient
params.show_iter_info = true;		% Flag for Showing Iteration Information

%% Calling PSO function
out = particleSwarm(problem, params);

pop = out.pop;
global_best = out.global_best;
best_costs = out.best_costs;
output = global_best.position;

%% Results
figure;
semilogy(output, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Costs');
