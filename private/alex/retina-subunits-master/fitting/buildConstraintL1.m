function out = buildConstraintL1(I_sc)

%-------------------------------------------
%
% out = buildConstraintL1(I_sc)
%
% build a matrix to implement an L1 constraint
%
% inputs
% I_sc -- matrix of binary assignments from
%         C stimuli to S subunits
%
% outputs
% out -- matrix that computes the sum of each
%        group of parameters corresponding to
%        one row (for use with fmincon)
% 
% freeman, 3-22-2012
%-------------------------------------------

%%

nSub = size(I_sc,1);
nPrs = length(find(I_sc));
I_sc_count = bsxfun(@times,I_sc,[1:nSub]');
foo = I_sc_count(find(I_sc));

counted = [1:nSub];
out = [];
for i=1:length(foo)
    if sum(counted==foo(i))
        newRow = zeros(1,nPrs);
        newRow(foo==foo(i)) = 1;
        out = [out; newRow];
        counted(foo(i)) = 0;
    end
end