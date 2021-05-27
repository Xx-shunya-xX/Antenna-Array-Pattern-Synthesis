%% Mutation Function Script
% r = true for real coded GA
% r = false for binary GA
function y = mutate(x, Mu, real, sigma)
	flag = (rand(size(x)) < Mu);
	y = x;
	if(real == true)
		r = randn(size(x));
		y(flag) = x(flag) + sigma*r(flag);
	else
		y(flag) = 1 - x(flag);
	end
end
