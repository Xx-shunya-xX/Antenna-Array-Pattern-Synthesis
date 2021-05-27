%% Using PSO with min max fitness function for linear antenna array
clc;
clear all;
close all;

%% Parameters for antenna array
f = 2.4e9;						% Frequency
c = 3e8;						% Speed of light
lambda = c / f;					% Wavelength
d = lambda / 2;					% distance between two array elements
k = (2 * pi) / lambda;			% Wave Number
theta = (0:0.1:180);			% 0 <= theta <= 76 && 104 <= theta <= 180
psi_n = 0;						% Phase Excitation

w_max = 0.9;
w_min = 0.5;

%% Problem Definition
problem.cost_func = @(x) minMax(x);		% Cost Function
problem.n_var = 5;						% Number of Unknown (Decision) Variables
problem.var_min = 0.1;				% Lower Bound of Decision
problem.var_max = 1;				% Upper Bound of Decision

%% Parameters of PSO
params.max_it = 100;				% Maximum Number of Iterations
params.n_pop = 40;					% Population (Swarm) Size
params.w = w_min - (w_max - w_min) / params.max_it;						% Inertia Coefficient
params.w_damp = 1;					% Damping Ration of Inertia Weight
params.c1 = 2;				% Personal Acceleration Coefficient
params.c2 = 2;				% Social Acceleration Coefficient
params.show_iter_info = true;		% Flag for Showing Iteration Information

%% Calling PSO function
out = particleSwarm(problem, params);
pop = out.pop;
global_best = out.global_best;
best_costs = out.best_costs;
I = global_best.position .'/ max(global_best.position);


%% Array Factor of antenna array
AF = arrayFactor(I, k, d, theta, psi_n)

AF_n = 20 * log10(AF/max(AF))

%% Results
plot(theta, AF_n);
xlabel('Azimuth angle (deg)');
ylabel('Array Factor');
legend('PSO');
grid on;
