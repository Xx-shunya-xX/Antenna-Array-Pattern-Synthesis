%% Sort Population
function pop = sortPopulation(pop)
	[~, so] = sort([pop.cost]);
	pop = pop(so);
end
