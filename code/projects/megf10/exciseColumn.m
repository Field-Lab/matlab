function m = exciseColumn(m)
% removes any columns that contain Inf or Nans from matrix
m(:, any(isnan(m),1)) = [];
m(:, any(isinf(m),1)) = [];
end
