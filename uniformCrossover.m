%% Uniform Crossover Function Script
% r = true for real coded GA
% r = false for binary GA
% Did'nt used gamma parameter
function [y1, y2] = uniformCrossover(x1, x2, real)
	%gamma = unifrnd(-gamma, 1+gamma, size(x1));
	if(real == true)
		alpha = rand(size(x1));
	else
		alpha = randi([0, 1], size(x1));
	end
	y1 = alpha .* x1 + (1 - alpha) .* x2;
	y2 = alpha .* x2 + (1 - alpha) .* x1;
end
