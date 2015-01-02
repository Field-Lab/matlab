function M = zerodimtoempty(M)
% Convert an empty array with nonzero dimensions to []
%
% Certain operations like repmat(1, 0, 2) can produce an empty array with
% nonzero dimensions (a 0-by-2 array in this example case).  This can
% confuse operations like +, which doesn't know how to add a 0-by-2 array
% to a 0-by-0 array.  So this function will drop a 0-by-2 array to [].

if isempty(M), M = []; end