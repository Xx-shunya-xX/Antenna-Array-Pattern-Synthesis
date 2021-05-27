% Sine Cosine Algorithm
function out = sincosineAlgorithm(problem, params)

	count=0;
	t11=1;                            % Number of Monte Carlo runs
	for tt=1:t11
	di=[25,28];                       % Direction of signals in degrees
	dii=abs(di(1)-di(2))/2;           % Difference between incoming signals
	d1=length(di);                    % Length of signals

	fobj = @problem.cost_func;
	Max_iteration = params.max_it;
	SearchAgents_no = params.n_pop;
	dim = problem.n_var; 
	var_min = problem.var_min;
	var_max = problem.var_max;
	lb = [repmat(var_min, 1, dim)];
	ub = [repmat(var_max, 1, dim)];

	%Initialize the set of random solutions
	for i=1:dim
		ub_i=ub(i);
		lb_i=lb(i)
		X(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
	end
	Destination_position=zeros(1,dim);
	Destination_fitness=inf;
	Convergence_curve=zeros(1,Max_iteration);
	Objective_values = zeros(1,size(X,1));
	% Calculate the fitness of the first set and find the best one
	for i=1:size(X,1)
		Objective_values(1,i)=feval(fobj,(X(i,:)));
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
		t=t+1;
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
			 
			% Check if solutions go outside the search spaceand bring them back
			Flag4ub=X(i,:)>ub;
			Flag4lb=X(i,:)<lb;
			X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
			
			% Calculate the objective values
			Objective_values(1,i)=feval(fobj,(X(i,:)));
			
			% Update the destination if there is a better solution
			if Objective_values(1,i)<Destination_fitness
				Destination_position=X(i,:);
				Destination_fitness=Objective_values(1,i);
			end
		end
		
		Convergence_curve(t)=Destination_fitness;
		format long
		[t Destination_position Destination_fitness]
	z11(tt,:)=sort(abs(Destination_position));
	y1(tt,1)=abs(z11(tt,1)-di(1));
	y2(tt,1)=abs(z11(tt,2)-di(2));
	y11(tt,1)=(z11(tt,1)-di(1)).^2;
	y22(tt,1)=(z11(tt,2)-di(2)).^2;
	if y1(tt)<dii && y2(tt)<dii
	   y0(tt,1)=1;
	   count=count+1;
	  else
		y0(tt,1)=0;
	end
	end

	end
	x1=sum(y11);
	x2=sum(y22);
	x3=x1+x2;
	format long;
	count;
	RMSE=sqrt(x3/(d1*t11));

	out.rmse = RMSE;
	out.global_best = Destination_position;
	out.best_cost = Destination_fitness;


end
