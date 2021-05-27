%% Using PSO with min max fitness function for linear antenna array
clc;
clear all;
close all;

%% Parameters for antenna array
f = input('Enter freq : ');							% Frequency
c = 3e8;											% Speed of light
lambda = c / f;										% Wavelength
d = lambda / 2;										% distance between two array elements
k = (2 * pi) / lambda;								% Wave Number
theta = (0:0.1:180);								% 0 <= theta <= 76 && 104 <= theta <= 180
N = input('Enter number of unknown variables : ');
pop = input('Enter number of population size : ');
psi_n = 0;											% Phase Excitation
max_it = input('Enter number of iterations : ');	% Maximum Number of Iterations
var_min = input('Enter lower bound : ');
var_max = input('Enter upper bound : ');

%% Problem Definition
problem.cost_func = @(x) minMax(x);		% Cost Function
problem.n_var = N;						% Number of Unknown (Decision) Variables
problem.var_min = var_min;				% Lower Bound of Decision
problem.var_max = var_max;				% Upper Bound of Decision

%% Parameters of PSO
% Using Clerc and Kennedy Constrictions
kappa = 1;							% 0 <= k <= 1
phi1 = 2.05;
phi2 = 2.05;
phi = phi1 + phi2;
chi = 2 * kappa / abs(2 - phi - sqrt(power(phi, 2) - 4 * phi));

params.max_it = max_it;				% Maximum Number of Iterations
params.n_pop = pop;					% Population (Swarm) Size
params.w = chi;						% Inertia Coefficient
params.w_damp = 1;					% Damping Ration of Inertia Weight
params.c1 = chi * phi1;				% Personal Acceleration Coefficient
params.c2 = chi * phi2;				% Social Acceleration Coefficient
params.show_iter_info = true;		% Flag for Showing Iteration Information

%% Calling PSO function
out = particleSwarm(problem, params);
pop = out.pop;
global_best = out.global_best;
best_costs = out.best_costs;
I = best_costs;

%% Array Factor of antenna array
for(n = 1 : length(theta))
	AF = 2 * I * cos(k * d * cos(theta(n)) + psi_n);
end

AF_n = 20 * log10(AF/max(AF));

%% Results
plot(theta, AF_n);
xlabel('Azimuth angle (deg)');
ylabel('Array Factor');
legend('PSO');
grid on;
