%% Particle Swarm Optimization Function Script
function out = particleSwarm(problem, params)
	%% Problem Definition
	cost_func = problem.cost_func;		% Cost Function
	n_var = problem.n_var;				% Number of Unknown (Decision) Variables
	var_size = [1 n_var];				% Matrix Size of Decision Variables
	var_min = problem.var_min;			% Lower Bound of Decision
	var_max = problem.var_max;			% Upper Bound of Decision

	%% Parameters of PSO
	max_it = params.max_it;				% Maximum Number of Iterations
	n_pop = params.n_pop;				% Population (Swarm) Size
	w = params.w;						% Inertia Coefficient
	w_damp = params.w_damp;				% Damping Ration of Inertia Weight
	c1 = params.c1;						% Personal Acceleration Coefficient
	c2 = params.c2;						% Social Acceleration Coefficient

	% Flag for Showing Iteration Information
	show_iter_info = params.show_iter_info;

	% New Addition to PSO to limit velocity to remain inside the defined limit
	max_velocity = 0.2 * (var_max - var_min);
	min_velocity = -max_velocity;

	%% Initialization
	% The Particle Template
	empty_particle.position = [];
	empty_particle.velocity = [];
	empty_particle.cost = [];
	empty_particle.best.position = [];
	empty_particle.best.cost = [];

	% Create Population Array
	particle = repmat(empty_particle, n_pop, 1);

	% Initialize Global Best
	global_best.cost = inf;

	% Initialize Population Members
	for(i = 1:n_pop)
		% Generate Random Solutions
		particle(i).position = unifrnd(var_min, var_max, var_size);

		% Initialize Velocity
		particle(i).velocity = zeros(var_size);

		% Evalution
		particle(i).cost = cost_func(particle(i).position);

		% Update the Personal Best
		particle(i).best.position = particle(i).position;
		particle(i).best.cost = particle(i).cost;

		% Update Global Best
		if(particle(i).best.cost < global_best.cost)
			global_best = particle(i).best;
		end
	end

	% Array to Hold Best Cost Values on Each Iteration
	best_costs = zeros(max_it, 1);

	%% Main Loop of PSO

	for(i = 1:max_it)
		for(j = 1:n_pop)
			% vij(t) = inertia term + cognitive component + social component
			% where; c1 & c2 are acceleration coefficients
			rand_number = rand(var_size);

			% Update Velocity
			particle(j).velocity = (w * particle(j).velocity) + (rand_number * c1 .* (particle(j).best.position - particle(j).position)) + (rand_number * c2 .* (global_best.position - particle(j).position));

			% Apply Velocity Limits
			particle(j).velocity = max(particle(j).velocity, min_velocity);
			particle(j).velocity = min(particle(j).velocity, max_velocity);

			% Update Position
			particle(j).position = particle(j).position + particle(j).velocity;

			% Apply upper bound and lower bound limits
			particle(j).position = max(particle(j).position, var_min);
			particle(j).position = min(particle(j).position, var_max);

			% Evaluation
			particle(j).cost = cost_func(particle(j).position);

			% Update the Personal Best
			if(particle(j).cost < particle(j).best.cost)
				particle(j).best.position = particle(j).position;
				particle(j).best.cost = particle(j).cost;

				% Update Global Best
				if(particle(j).best.cost < global_best.cost)
					global_best = particle(j).best;
				end
			end
		end
		% Store the Best Cost Value
		best_costs(i) = global_best.cost;

		%Display Iteration Information
		if(show_iter_info)
			disp(['Iteration ' num2str(i) ' : Best Cost = ' num2str(best_costs(i))]);
			disp('Global Best : ');
			disp(global_best);
		end

		% Damping Inertia Coefficient
		w = w * w_damp;
	end

	out.pop = particle;
	out.global_best = global_best;
	out.best_costs = best_costs;
end
