function cost = cone_spacing_cost_function(centers)
% For simulated annealing

distances = ipdm(centers);
costs = (1 ./ distances).^2;
costs(costs == Inf) = 0;
cost = sum(costs(:)) ./ 2;