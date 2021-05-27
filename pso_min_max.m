% pso demo with constriction coefficients
clc;
clear all;
close all;

array_result = zeros(15, 5);

%% Problem Definition
problem.cost_func = @(x) minMax(x);	% Cost Function
problem.n_var = 5;					% Number of Unknown (Decision) Variables
problem.var_min = 0.00001;				% Lower Bound of Decision
problem.var_max = 0.99999;				% Upper Bound of Decision

%% Parameters of PSO
% Clerk and kennedy constrictions
kappa = 1;							% 0 <= k <= 1
phi1 = 2.05;
phi2 = 2.05;
phi = phi1 + phi2;
chi = 2 * kappa / abs(2 - phi - sqrt(power(phi, 2) - 4 * phi));

w_max = 0.9;
w_min = 0.4;
params.max_it = 1000;					% Maximum Number of Iterations
params.n_pop = 40;					% Population (Swarm) Size
params.w = w_min - (w_max - w_min) / params.max_it;						% Inertia Coefficient
params.w_damp = 1;					% Damping Ration of Inertia Weight
params.c1 = chi * phi1;				% Personal Acceleration Coefficient
params.c2 = chi * phi2;				% Social Acceleration Coefficient
%params.c1 = 2;
%params.c2 = 2;
params.show_iter_info = true;		% Flag for Showing Iteration Information

%% Calling PSO function
for(itr = 1 : 15)
	out = particleSwarm(problem, params);
	pop = out.pop;
	global_best = out.global_best;
	best_costs = out.best_costs;
	sorted_result = sort(global_best.position ./ max(global_best.position), 'descend');
	array_result(itr, :) = sorted_result;
end

result = median(array_result);

%% Results
disp('Global Best positions : ');
disp(result);
figure;
plot(result, 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Costs');
