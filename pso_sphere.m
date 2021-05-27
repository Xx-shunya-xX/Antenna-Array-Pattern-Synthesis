clc;
clear all;
close all;

%% Problem Definition
problem.cost_func = @(x) Sphere(x);	% Cost Function
problem.n_var = 2;					% Number of Unknown (Decision) Variables
problem.var_min = -3;				% Lower Bound of Decision
problem.var_max = 3;				% Upper Bound of Decision

%% Parameters of PSO
params.max_it = 100;				% Maximum Number of Iterations
params.n_pop = 50;					% Population (Swarm) Size
params.w = 1;						% Inertia Coefficient
params.w_damp = 0.99;				% Damping Ration of Inertia Weight
params.c1 = 2;						% Personal Acceleration Coefficient
params.c2 = 2;						% Social Acceleration Coefficient
params.show_iter_info = true;		% Flag for Showing Iteration Information

%% Calling PSO function
out = particleSwarm(problem, params);

pop = out.pop;
global_best = out.global_best;
best_costs = out.best_costs;

%% Results
figure;
semilogy(best_costs, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Costs');
