function [y1, y2] = randomCrossover(x1, x2)
	m = randi([1, 3]);
	switch m
		case 1
			[y1, y2] = singlePointCrossover(x1, x2);

		case 2
			[y1, y2] = doublePointCrossover(x1, x2);

		otherwise
			[y1, y2] = uniformCrossover(x1, x2, false);
	end
end
