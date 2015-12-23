function norm_sum = sum_norms(response_subspace, init_weights)

% an arbitrary set of weights

subspace_size = size(response_subspace,2);

% form linear combinations of the response space vectors to get a new basis
% set


linear_combinations = response_subspace * init_weights;


l_one_norms = ones(subspace_size,1);
for dm = 1:size(linear_combinations,2);
    l_one_norms(dm) = norm(linear_combinations(:,dm),1);
end

norm_sum = prod(l_one_norms);








