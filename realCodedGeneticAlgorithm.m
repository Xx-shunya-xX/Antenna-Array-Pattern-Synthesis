%% Real Coded Genetic Algorithm Function Script
% Did'nt used gamma parameter for unifromCrossover for unifrnd initialization
function out = realCodedGeneticAlgorithm(problem, params)
	% Problem
	cost_func = problem.cost_func;
	n_var = problem.n_var;
	var_size = [1, n_var];
	var_min = problem.var_min;
	var_max = problem.var_max;

	% Params
	max_it = params.max_it;
	n_pop = params.n_pop;
	beta = params.beta;									% Beta value for probs
	percent_child = params.percent_child;				% Percentage of Children
	Mu = params.Mu;
	numb_child = round(percent_child * n_pop / 2) * 2;	% Number of Children
	sigma = params.sigma;								% Sigma value for mutate
	show_iter_info = params.show_iter_info;

	% Template for Empty Individuals
	empty_individual.position = [];
	empty_individual.cost = [];

	% Best Solution Ever Found
	global_best.cost = inf;

	% Initialization
	pop = repmat(empty_individual, n_pop, 1);
	for(i = 1:n_pop)

		% Generate Uniform Random Solution
		pop(i).position = unifrnd(var_min, var_max, var_size);

		% Evaluate Solution
		pop(i).cost = cost_func(pop(i).position);

		% Compare Solution to Best Solution Ever Found
		if(pop(i).cost < global_best.cost)
			global_best = pop(i);
		end
	end
	
	% Best Cost of Iterations
	best_cost = nan(max_it, 1);

	% Main Loop
	for(i = 1:max_it)

		% Selection Probabilities
		c = [pop.cost];
		avgc = mean(c);
		if(avgc ~= 0)
			c = c / avgc;
		end
		probs = exp(-beta * c);

		% Initialize Offsprings Population
		pop_c = repmat(empty_individual, numb_child / 2, 1);

		% Crossover
		for(j = 1:numb_child / 2)
			% Select Parents
			p1 = pop(rouletteWheelSelection(probs));
			p2 = pop(rouletteWheelSelection(probs));

			% Perform Crossover
			[pop_c(j, 1).position, pop_c(j, 2).position] = uniformCrossover(p1.position, p2.position, true);
		end

		% Convert pop_c to Single-Column Matrix
		pop_c = pop_c(:);

		% Mutation
		for(k = 1 : numb_child)
			% Perform Mutation
			pop_c(k).position = mutate(pop_c(k).position, Mu, true, sigma);

			% Upper and Lower Bounds
			pop_c(k).position = max(pop_c(k).position, var_min);
			pop_c(k).position = min(pop_c(k).position, var_max);

			% Evaluate
			pop_c(k).cost = cost_func(pop_c(k).position);

			% Compare Solution to Best Solution Ever Found
			if(pop_c(k).cost < global_best.cost)
				global_best = pop_c(k);
			end
		end

		% Merge and Sort Populations
		pop = sortPopulation([pop; pop_c]);

		% Remove Extra Individuals
		pop = pop(1 : n_pop);

		% Update Best Cost of Iteration
		best_cost(i) = global_best.cost;

		% Display Iteration Information
		if(show_iter_info)
			disp(['Iteration ' num2str(i) ' : Best Cost = ' num2str(best_cost(i))]);
			disp('Global Best : ');
			disp(global_best);
		end
	end

	% Return Result
	out.pop = pop;
	out.global_best = global_best;
	out.best_cost = best_cost;
end
