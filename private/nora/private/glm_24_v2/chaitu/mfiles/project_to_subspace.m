% Project a vector v to span of columns of B with lambda L2-regularizer
function Pv = project_to_subspace(B,v,lambda)

if (nargin < 3)
    lambda = 0;
end


Pv = (B'*B + lambda.*eye(size(B,2)))\(B'*v);