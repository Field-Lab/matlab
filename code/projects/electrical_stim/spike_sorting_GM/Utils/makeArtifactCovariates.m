function [X Xj] = makeArtifactCovariates(T,J,I,varargin)
X = sparse(0,0);
if(nargin == 3)
    E = 1;
else
    E = varargin{4};
end

for j = 1:J
       ind    = sparse(1,J);
       ind(j) = 1;
       Xj{j}  = kron(ind,speye(E*T));
       X      = [X;repmat(Xj{j},I(j),1)];
end