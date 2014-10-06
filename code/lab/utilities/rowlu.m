function results = rowlu(data, lookups)
% ROWLU    Get array elements according to coordinates given as row vectors
% example: d = [1 2 3; 4 5 6; 7 8 9];
%          lu = [1 1; 2 2; 3 3];
%          rowlu(d, lu)
%          >> [1 5 9]

% ALGORITHM STRATEGY IS TO CONVERT EACH LOOKUP ROW TO A 1D INDEX

if size(lookups, 2) ~= numel(size(data))
    error('LOOKUPS must have as many columns as the dimensionality of DATA.')
end

multipliers = cumprod(size(data));
multipliers = [1 multipliers(1:end-1)];
multmat = repmat(multipliers, size(lookups, 1), 1);

zb_lookups = lookups - 1; % zero-base the lookups
indices = sum(multmat .* zb_lookups, 2);
indices = indices + 1; % return to one-base

results = data(indices);