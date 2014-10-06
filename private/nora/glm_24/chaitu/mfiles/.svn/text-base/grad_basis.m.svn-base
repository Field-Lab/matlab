% Convolve the spike trains with a set of basis functions

function A = grad_basis(D,basis,offset)

if (nargin < 3)
    offset = 0;
end

[T ntrains] = size(D);
%[m nbasis] = size(basis);

A = cell(ntrains,1);

for j=1:ntrains
    A{j} = spike_convmex(T,find(D(:,j)),basis)';
    A{j} = A{j}(:,offset+1:end);
end        


