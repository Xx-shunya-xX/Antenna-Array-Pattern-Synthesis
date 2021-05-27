%% Using PSO with min max fitness function for linear antenna array
clc;
clear all;
close all;

%% Parameters for antenna array
f = input('Enter freq : ');					% Frequency
c = 3e8;									% Speed of light
lambda = c / f;								% Wavelength
d = lambda / 2;								% distance between two array elements
k = (2 * pi) / lambda;						% Wave Number
theta_1 = (0 : 0.1 : 76);
theta_2 = (104 : 0.1 : 180);
theta = [theta_1, theta_2];					% 0 <= theta <= 76 && 104 <= theta <= 180
N = input('Enter number of elements : ');
psi = 0;									% Phase Excitation
max_it = 5;									% Maximum Number of Iterations
I1 = 1;
I1 = repmat(I1, 1, 10);

AF = arrayFactor(N, I1, k, d, theta, psi);

%% Problem Definition
problem.cost_func = @(x) minMax(x);		% Cost Function
problem.n_var = N;						% Number of Unknown (Decision) Variables
problem.var_min = 0;					% Lower Bound of Decision
problem.var_max = 1;					% Upper Bound of Decision

%% Parameters of PSO
% Using Clerc and Kennedy Constrictions
kappa = 1;							% 0 <= k <= 1
phi1 = 2.05;
phi2 = 2.05;
phi = phi1 + phi2;
chi = 2 * kappa / abs(2 - phi - sqrt(power(phi, 2) - 4 * phi));

params.max_it = N;					% Maximum Number of Iterations
params.n_pop = N;					% Population (Swarm) Size
params.w = chi;						% Inertia Coefficient
params.w_damp = 1;					% Damping Ration of Inertia Weight
params.c1 = chi * phi1;				% Personal Acceleration Coefficient
params.c2 = chi * phi2;				% Social Acceleration Coefficient
params.show_iter_info = false;		% Flag for Showing Iteration Information

%% Calling PSO function
out = particleSwarm(problem, params);
pop = out.pop;
global_best = out.global_best;
bests_costs = out.best_costs;
I2 = out.best_costs;

AF2 = arrayFactor(N, I2, k, d, theta, psi);

%% Results
plot(theta, AF, theta, AF2);
xlabel('Azimuth angle (deg)');
ylabel('Array Factor');
legend('AF', 'PSO optimized AF');
grid on;
