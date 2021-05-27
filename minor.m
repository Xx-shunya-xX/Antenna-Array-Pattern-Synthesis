clc;
clear all;
close all;

format long;

%% Set Metaheuristic Algorithms ON/OFF
optimize_results = true;			% Optmize the results or not, if false, show only convenction and proposed results
optimize_with_pso = true;			% Optimize with Particle Swarm Optimization
optimize_with_rga = true;			% Optimize with Real Coded Genetical Algorithm
optimize_with_tsa = true;			% Optimize with Tunicate Swarm Algorithm
optimize_with_sca = true;			% Optimize with Sine Cosine Algorithm

separated_graphs = true;			% If you want separated graphs for less visible confusion
params.show_iter_info = true;		% Flag for Showing Iteration Information
save_results = false;

%% Set Axis for graph plots
axis_x_from = 0;
axis_x_to = 180;
axis_y_from = -40;
axis_y_to = 0;

last_fig_number = 1;	% Variable to store figure counts

%% Antenna Parameters
f = 2.4e9;				% Frequency (Default : 2.4 GHz)
c = 3e8;				% Speed of light
lambda = c / f;			% Wavelength
k = pi / 2;				% Wave Number
psi = 0;				% Amplitude Excitation
d = lambda / 2;			% Distance between two elements of array
theta = (0:0.1:180);	% Angle Division
w_max = 0.9;			% Weight Max
w_min = 0.5;			% Weight Min
N = 5;					% Number of variables
no_of_runs = 15;		% Number of runs

% Clerk and kennedy constrictions
kappa = 1;							% 0 <= k <= 1
phi1 = 2.05;
phi2 = 2.05;
phi = phi1 + phi2;
chi = 2 * kappa / abs(2 - phi - sqrt(power(phi, 2) - 4 * phi));

%% Problem Definition
problem.cost_func = @(x) minMax(x);		% Cost Function
problem.n_var = N;						% Number of Unknown (Decision) Variables
problem.dimensions = problem.n_var;
problem.var_min = 0.000001;				% Lower Bound of Decision
problem.var_max = 0.999999;				% Upper Bound of Decision

%% Parameters of PSO
params.max_it = 1000;										% Maximum Number of Iterations
params.n_pop = 40;										% Population (Swarm) Size
params.w = w_min - (w_max - w_min) / params.max_it;		% Inertia Coefficient
params.w_damp = 1;										% Damping Ration of Inertia Weight
params.c1 = chi * phi1;											% Personal Acceleration Coefficient
params.c2 = chi * phi2;											% Social Acceleration Coefficient
params.percent_child = 1;
params.Mu = 0.00001;
params.beta = 1;
params.sigma = 0.2;

array_af_pso = zeros(no_of_runs, N);
array_af_rga = zeros(no_of_runs, N);
array_af_tsa = zeros(no_of_runs, N);
array_af_sca = zeros(no_of_runs, N);

%% Conventional & Proposed PSO optimized result
conven_amp_excit = ones(5, 1);
proposed_amp_excit = [1; 0.9018; 0.72759; 0.51502; 0.4159];

%% Calling all the optimization techniques for number of iterations
if(optimize_results == false)
	optimize_with_pso = false;			% Optimize with Particle Swarm Optimization
	optimize_with_rga = false;			% Optimize with Real Coded Genetical Algorithm
	optimize_with_tsa = false;			% Optimize with Tunicate Swarm Algorithm
	optimize_with_sca = false;			% Optimize with Sine Cosine Algorithm
end

for(nor = 1 : no_of_runs)
	if(optimize_with_pso)
		% Calling PSO function (AF_pso)
		disp('----------------------------> Optimizing Amplitude excitation by PSO');
		disp(['Number of global iterations running on : ' num2str(nor)])
		out_pso = particleSwarm(problem, params);
		I_pso = sort(out_pso.global_best.position .'/ max(out_pso.global_best.position), 'descend');

		array_af_pso(nor, :) = I_pso;
	end

	if(optimize_with_rga)
		% Calling Real Genetic Algorithm (AF_rga)
		disp('----------------------------> Optimizing Amplitude excitation by RGA');
		disp(['Number of global iterations running on : ' num2str(nor)])
		out_rga = realCodedGeneticAlgorithm(problem, params);
		I_rga = sort(out_rga.global_best.position .'/ max(out_rga.global_best.position), 'descend');

		array_af_rga(nor, :) = I_rga;
	end

	if(optimize_with_tsa)
		% TSA algo
		disp('----------------------------> Optimizing Amplitude excitation by TSA');
		disp(['Number of global iterations running on : ' num2str(nor)])
		out_tsa = tunicateAlgorithm(problem, params);
		I_tsa = sort(out_tsa.global_best .'/ max(out_tsa.global_best), 'descend');

		array_af_tsa(nor, :) = I_tsa;
	end

	if(optimize_with_sca)
		% SCA algo
		disp('----------------------------> Optimizing Amplitude excitation by SCA');
		disp(['Number of global iterations running on : ' num2str(nor)])
		out_sca = sineCosineAlgorithm(problem, params);
		I_sca = sort(out_sca.global_best .'/ max(out_sca.global_best), 'descend');

		array_af_sca(nor, :)= I_sca;
	end
end % No of runs

%% Calculating Median & Standard Deviation of optimized results
if(optimize_with_pso)
	I_pso_median = median(array_af_pso);
	I_pso_std = std(array_af_pso);
end

if(optimize_with_rga)
	I_rga_median = median(array_af_rga);
	I_rga_std = std(array_af_rga);
end

if(optimize_with_tsa)
	I_tsa_median = median(array_af_tsa);
	I_tsa_std = std(array_af_tsa);
end

if(optimize_with_sca)
	I_sca_median = median(array_af_sca);
	I_sca_std = std(array_af_sca);
end

%% Calculating Array Factor
AF_conven = arrayFactor(conven_amp_excit, k, 5, theta, psi);
AF_proposed = arrayFactor(proposed_amp_excit, k, 5, theta, psi);

if(optimize_with_pso)
	AF_pso = arrayFactor(I_pso_median.', k, N, theta, psi);
end
if(optimize_with_rga)
	AF_rga = arrayFactor(I_rga_median.', k, N, theta, psi);
end
if(optimize_with_tsa)
	AF_tsa = arrayFactor(I_tsa_median.', k, N, theta, psi);
end
if(optimize_with_sca)
	AF_sca = arrayFactor(I_sca_median.', k, N, theta, psi);
end

%% Displaying Results
disp('Conventional values : ');
disp(conven_amp_excit);
disp('Proposed values : ');
disp(proposed_amp_excit);

if(optimize_with_pso)
	disp('PSO : ');
	disp(out_pso.global_best.position);
	disp('Amplitude Excitation by PSO : ');
	disp(I_pso);
end

if(optimize_with_rga)
	disp('Real Coded GA : ');
	disp(out_rga.global_best.position);
	disp('Amplitude Excitation by Real Coded GA : ');
	disp(I_rga);
end

if(optimize_with_tsa)
	disp('Global best by TSA : ');
	disp(out_tsa.global_best);
	disp('Amplitude Excitation by TSA : ');
	disp(I_tsa);
end

if(optimize_with_sca)
	disp('Global best by SCA : ');
	disp(out_sca.global_best);
	disp('Amplitude Excitation by SCA : ');
	disp(I_sca);
end

%% Plotting Graph
figure(last_fig_number)
last_fig_number = last_fig_number + 1;
plot(theta, 20 * log10(AF_conven / max(AF_conven)), 'k-.', 'DisplayName', 'Conv')
hold on
plot(theta, 20 * log10(AF_proposed / max(AF_proposed)), 'Marker', '.', 'MarkerFaceColor', [0.4940 0.1840 0.5560], 'DisplayName', 'Proposed PSO')
hold on

if(optimize_with_pso)
	plot(theta, 20 * log10(AF_pso / max(AF_pso)), 'r', 'DisplayName', 'Mine PSO')
	hold on
end

if(optimize_with_rga)
	plot(theta, 20 * log10(AF_rga / max(AF_rga)), 'g', 'DisplayName', 'Real GA')
	hold on
end

if(optimize_with_tsa)
	plot(theta, 20 * log10(AF_tsa / max(AF_tsa)), 'b', 'DisplayName', 'TSA')
	hold on
end

if(optimize_with_sca)
	plot(theta, 20 * log10(AF_sca / max(AF_sca)), 'm', 'DisplayName', 'SCA')
	hold off
end

axis([axis_x_from axis_x_to axis_y_from axis_y_to])
ylabel('Gain (dB)')
xlabel('Azimuth angle (deg)')
legend
grid on

%% Element Graph
figure(last_fig_number)
last_fig_number = last_fig_number + 1;

elements_neg = (-5 : -1);
elements = [elements_neg, fliplr(abs(elements_neg))];

%% Calculating Normalized amplitudes
rr_conven = conven_amp_excit;
norm_amp_conven = [rr_conven(5) rr_conven(4) rr_conven(3) rr_conven(2) rr_conven(1) rr_conven(1) rr_conven(2) rr_conven(3) rr_conven(4) rr_conven(5)];

rr_proposed = proposed_amp_excit;
norm_amp_proposed = [rr_proposed(5) rr_proposed(4) rr_proposed(3) rr_proposed(2) rr_proposed(1) rr_proposed(1) rr_proposed(2) rr_proposed(3) rr_proposed(4) rr_proposed(5)];

if(optimize_with_pso)
	% PSO
	rr_pso = I_pso;
	norm_amp_pso = [rr_pso(5) rr_pso(4) rr_pso(3) rr_pso(2) rr_pso(1) rr_pso(1) rr_pso(2) rr_pso(3) rr_pso(4) rr_pso(5)];
end

if(optimize_with_rga)
	% RGA
	rr_rga = I_rga;
	norm_amp_rga = [rr_rga(5) rr_rga(4) rr_rga(3) rr_rga(2) rr_rga(1) rr_rga(1) rr_rga(2) rr_rga(3) rr_rga(4) rr_rga(5)];
end

if(optimize_with_tsa)
	% TSA
	rr_tsa = I_tsa;
	norm_amp_tsa = [rr_tsa(5) rr_tsa(4) rr_tsa(3) rr_tsa(2) rr_tsa(1) rr_tsa(1) rr_tsa(2) rr_tsa(3) rr_tsa(4) rr_tsa(5)];
end

if(optimize_with_sca)
	% SCA
	rr_sca = I_sca;
	norm_amp_sca = [rr_sca(5) rr_sca(4) rr_sca(3) rr_sca(2) rr_sca(1) rr_sca(1) rr_sca(2) rr_sca(3) rr_sca(4) rr_sca(5)];
end

% Plotting graph
plot(elements, norm_amp_conven, 'k-.', 'DisplayName', 'Conv');
hold on

plot(elements, norm_amp_proposed, 'k--', 'DisplayName', 'Proposed');
hold on

if(optimize_with_pso)
	plot(elements, norm_amp_pso, 'r-o', 'DisplayName', 'Mine PSO');
	hold on
end

if(optimize_with_rga)
	plot(elements, norm_amp_rga, 'g-o', 'DisplayName', 'Real GA');
	hold on
end

if(optimize_with_tsa)
	plot(elements, norm_amp_tsa, 'b-o', 'DisplayName', 'TSA');
	hold on
end

if(optimize_with_sca)
	plot(elements, norm_amp_sca, 'm-o', 'DisplayName', 'SCA');
	hold off
end

xlabel('Elements #')
ylabel('Normalized amplitude')
axis([-5 5 0.4 1.1])
legend
title('Amplitude vs Elements Graph')

% ----------------------------- >>>> Separated Elements graphs
if(separated_graphs)
end

% Separated AF graph
if(separated_graphs)
	figure(last_fig_number)
	last_fig_number = last_fig_number + 1;

	if(optimize_with_pso)
		subplot(2,2,1)
		plot(theta, 20 * log10(AF_conven / max(AF_conven)), 'k-.', 'DisplayName', 'Conv')
		hold on
		plot(theta, 20 * log10(AF_proposed / max(AF_proposed)), 'Marker', '.', 'MarkerFaceColor', [0.4940 0.1840 0.5560], 'DisplayName', 'Proposed PSO')
		hold on
		plot(theta, 20 * log10(AF_pso / max(AF_pso)), 'r', 'DisplayName', 'PSO')
		hold off
		axis([axis_x_from axis_x_to axis_y_from axis_y_to])
		ylabel('Gain (dB)')
		xlabel('Azimuth angle (deg)')
		legend('Conv', 'Proposed', 'PSO')
		grid on
	end

	if(optimize_with_rga)
		subplot(2,2,2)
		plot(theta, 20 * log10(AF_conven / max(AF_conven)), 'k-.', 'DisplayName', 'Conv')
		hold on
		plot(theta, 20 * log10(AF_proposed / max(AF_proposed)), 'Marker', '.', 'MarkerFaceColor', [0.4940 0.1840 0.5560], 'DisplayName', 'Proposed PSO')
		hold on
		plot(theta, 20 * log10(AF_rga / max(AF_rga)), 'g', 'DisplayName', 'RGA')
		hold off
		axis([axis_x_from axis_x_to axis_y_from axis_y_to])
		ylabel('Gain (dB)')
		xlabel('Azimuth angle (deg)')
		legend('Conv', 'Proposed', 'RGA')
		grid on
	end

	if(optimize_with_tsa)
		subplot(2,2,3)
		plot(theta, 20 * log10(AF_conven / max(AF_conven)), 'k-.', 'DisplayName', 'Conv')
		hold on
		plot(theta, 20 * log10(AF_proposed / max(AF_proposed)), 'Marker', '.', 'MarkerFaceColor', [0.4940 0.1840 0.5560], 'DisplayName', 'Proposed PSO')
		hold on
		plot(theta, 20 * log10(AF_tsa / max(AF_tsa)), 'b', 'DisplayName', 'TSA')
		hold off
		axis([axis_x_from axis_x_to axis_y_from axis_y_to])
		ylabel('Gain (dB)')
		xlabel('Azimuth angle (deg)')
		legend('Conv', 'Proposed', 'TSA')
		grid on
	end

	if(optimize_with_sca)
		subplot(2,2,4)
		plot(theta, 20 * log10(AF_conven / max(AF_conven)), 'k-.', 'DisplayName', 'Conv')
		hold on
		plot(theta, 20 * log10(AF_proposed / max(AF_proposed)), 'Marker', '.', 'MarkerFaceColor', [0.4940 0.1840 0.5560], 'DisplayName', 'Proposed PSO')
		hold on
		plot(theta, 20 * log10(AF_sca / max(AF_sca)), 'm', 'DisplayName', 'SCA')
		hold off
		axis([axis_x_from axis_x_to axis_y_from axis_y_to])
		ylabel('Gain (dB)')
		xlabel('Azimuth angle (deg)')
		legend('Conv', 'Proposed', 'SCA')
		grid on
	end
end

%% Save results
if(save_results)
	result_path = "D:/Google Drive/College/SEM 6/minor-II/END/";
	result_counter_path = result_path + "result_counter";
	rc = load(result_counter_path);
	res_count = rc.res_count + 1;
	result_folder = ['results' num2str(res_count)];
	result_path = result_path + result_folder;
	delete(result_counter_path + ".mat");
	save(result_counter_path, "res_count");
	mkdir(result_path);
	graph_path = result_path + "/";
	%saveas(graph_power_pattern, graph_path + "graph_power_pattern", 'png');
	%saveas(graph_elements, graph_path + "graph_elements", 'png');
	%saveas(graph_separated, graph_path + "graph_separated", 'png');
end
