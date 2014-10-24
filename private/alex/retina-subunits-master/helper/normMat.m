function [out norms] = normMat(mat,type,dim)

%-------------------------------------------
%
% out = normMat(mat,type,dim)
%
% set the columns of a matrix to all have a norm of 1
%
% inputs:
% mat -- matrix that needs to have its columns normalized
% type -- type of norm (either 'L1' or 'L2')
% dim -- dimension along which to normalize
%
% freeman, 3-17-2012
%-------------------------------------------

% fixed a bug here 4-21-2031, badinds was in the wrong index, may have
% never had an effect if the number of rows was less than the number of
% columns...

if nargin < 2
    type = 'L2';
    dim = 1;
end

if nargin < 3
    dim = 1;
end

% get the norms
norms = getNorm(mat,type,dim);
badinds = norms==0;

% do fast point-wise division
out = mat;
if size(out,1) == 1 || size(out,2) == 1
    out = out / norms;
else
	if dim == 1
    out(:,~badinds) = bsxfun(@rdivide,mat(:,~badinds),norms(~badinds));
  elseif dim == 2
  	out(~badinds,:) = bsxfun(@rdivide,mat(~badinds,:),norms(~badinds));
  end	
end