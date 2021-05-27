% Tunicate Swarm Algorithm
function out = tunicateAlgorithm(problem, params)
tic; 
objective = @problem.cost_func;
Search_Agents = problem.n_var;
dimensions = Search_Agents;
Max_iterations = params.max_it;
show_iter_info = params.show_iter_info;
Lowerbound = [repmat(problem.var_min, 1, problem.dimensions)];
Upperbound = [repmat(problem.var_max, 1, problem.dimensions)];
Score = inf; 
for i = 1:dimensions
     ub_i = Upperbound(i);
     lb_i = Lowerbound(i);
     Positions(:,i) = rand(Search_Agents,1).*(ub_i-lb_i)+lb_i;
end
Convergence=zeros(1,Max_iterations);
t=0;

while t<Max_iterations
    for i=1:size(Positions,1)
        
    
        Flag4Upperbound=Positions(i,:)>Upperbound;
        Flag4Lowerbound=Positions(i,:)<Lowerbound;
        Positions(i,:)=(Positions(i,:).*(~(Flag4Upperbound+Flag4Lowerbound)))+Upperbound.*Flag4Upperbound+Lowerbound.*Flag4Lowerbound;
    
        fitness=objective(Positions(i,:));
    
		if fitness<Score 
            Score=fitness; 
            Position=Positions(i,:);
        end
        
    end
    
            xmin=1;
            xmax=4;
            xr=xmin+rand()*(xmax-xmin);
            xr=fix(xr);
  
   for i=1:size(Positions,1)
        for j=1:size(Positions,2)

           
             A1=((rand()+rand())-(2*rand()))/xr; 
                
            c2=rand();
        if(i==1)
        c3=rand(); 
        if(c3>=0)
            d_pos=abs(Position(j)-c2*Positions(i,j));
            Positions(i,j)=Position(j)+A1*d_pos;
        else
            d_pos=abs(Position(j)-c2*Positions(i,j));
        Positions(i,j)=Position(j)-A1*d_pos;
        
        
        
        end
        else
            
            
            c3=rand(); 
        if(c3>=0)
            d_pos=abs(Position(j)-c2*Positions(i,j));
            Pos(i,j)=Position(j)+A1*d_pos;
        else                          
        Pos(i,j)=Position(j)-A1*d_pos;
        
        
        end
             Positions(i,j)=(Pos(i,j)+Positions(i-1,j))/2;

			 if(show_iter_info)
				disp(['Iteration ' num2str(i)]);
				disp('Global Best Position : ');
				disp(Position);
			 end
        end
       
		end
            
        end
    
    t=t+1;
    Convergence(t)=Score;
	out.iterations = t;
	out.global_best = Position;
	out.best_cost = Score;
end
end
