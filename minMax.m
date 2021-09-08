function x = minMax(array_factor_phi)
	x = 0;
	theta_1 = (0 : 0.1 : 76);
	theta_2 = (104 : 0.1 : 180);
	theta = [theta_1, theta_2];		% 0 <= theta <= 76 && 104 <= theta <= 180
	L = (1 : length(theta));
	for(n = 1 : length(array_factor_phi))
		cost_func = 2 * array_factor_phi(n) * cos(((n * 2) - 1) * (pi/2) * cos((pi/180) * theta(L)));
		x = x + cost_func;
	end
	x = min(max(20 * log(abs(x))));
end
