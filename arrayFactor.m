%% Array Factor of antenna array
function AF = arrayFactor(amp_excit, wave_n, dist, theta, phase_excit)
	for(n = 1 : length(theta))
		AF(n)=abs((2 * (cos(wave_n * [1:2:((dist * 2) - 1)] * cos(theta(n) * (pi/180))))) * amp_excit);
	end
end
