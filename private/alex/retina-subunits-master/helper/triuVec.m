function out = triuVec(mat,k)

% get the upper triangular portion of a matrix using the triu function

if nargin < 2
    k = 0;
end

inds = logical(triu(ones(size(mat)),k));
out = mat(inds);