function out = sineCosineAlgorithm(problem, params)

	lb = problem.var_min;
	ub = problem.var_max;
	dim = problem.n_var;
	fobj = problem.cost_func;
	N = params.n_pop;
	Max_iteration = params.max_it;
	show_iter_info = params.show_iter_info;

	%Initialize the set of random solutions
	X=sineCosineInitialization(N,dim,ub,lb);

	Destination_position=zeros(1,dim);
	Destination_fitness=inf;

	Convergence_curve=zeros(1,Max_iteration);
	Objective_values = zeros(1,size(X,1));

	% Calculate the fitness of the first set and find the best one
	for i=1:size(X,1)
		Objective_values(1,i)=fobj(X(i,:));
		if i==1
			Destination_position=X(i,:);
			Destination_fitness=Objective_values(1,i);
		elseif Objective_values(1,i)<Destination_fitness
			Destination_position=X(i,:);
			Destination_fitness=Objective_values(1,i);
		end
		
		All_objective_values(1,i)=Objective_values(1,i);
	end

	%Main loop
	t=2; % start from the second iteration since the first iteration was dedicated to calculating the fitness
	while t<=Max_iteration
		
		% Eq. (3.4)
		a = 2;
		Max_iteration = Max_iteration;
		r1=a-t*((a)/Max_iteration); % r1 decreases linearly from a to 0
		
		% Update the position of solutions with respect to destination
		for i=1:size(X,1) % in i-th solution
			for j=1:size(X,2) % in j-th dimension
				
				% Update r2, r3, and r4 for Eq. (3.3)
				r2=(2*pi)*rand();
				r3=2*rand;
				r4=rand();
				
				% Eq. (3.3)
				if r4<0.5
					% Eq. (3.1)
					X(i,j)= X(i,j)+(r1*sin(r2)*abs(r3*Destination_position(j)-X(i,j)));
				else
					% Eq. (3.2)
					X(i,j)= X(i,j)+(r1*cos(r2)*abs(r3*Destination_position(j)-X(i,j)));
				end
				
			end
		end
		
		for i=1:size(X,1)
			 
			% Check if solutions go outside the search space and bring them back
			Flag4ub=X(i,:)>ub;
			Flag4lb=X(i,:)<lb;
			X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
			
			% Calculate the objective values
			Objective_values(1,i)=fobj(X(i,:));
			
			% Update the destination if there is a better solution
			if Objective_values(1,i)<Destination_fitness
				Destination_position=X(i,:);
				Destination_fitness=Objective_values(1,i);
			end
		end
		
		Convergence_curve(t)=Destination_fitness;
		
		% Display the iteration and best optimum obtained so far
		if(show_iter_info)
			display(['At iteration ', num2str(t), ' the optimum is ', num2str(Destination_fitness)]);
		end
		
		% Increase the iteration counter
		t=t+1;
	end
	out.global_best = Destination_position;
	out.best_cost = Destination_fitness;
end
