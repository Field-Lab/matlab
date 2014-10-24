function out = getNorm(mat,type,dim)

%-------------------------------------------
%
% out = getNorms(mat,type,dim)
%
% get the norm of a vector (or norms of a matrix's columns)
%
% inputs:
% mat - matrix (or vector)
% type - type of norm (either 'L1' or 'L2')
% dim - dimension to compute norms along
%
% freeman, 3-17-2012
%-------------------------------------------

if nargin < 2
    type = 'L2';
    dim = 1;
end

if nargin < 3
    dim = 1;
end

% compute the norms
switch type
    case 'L2'
        out = sqrt(sum(mat.^2,dim));
    case 'L1'
        out = (sum(abs(mat),dim));        
end