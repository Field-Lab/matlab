function simpleprint(M)
% SIMPLEPRINT   Simple print out of array M
%
% 2012-09 phli
%

for i = 1:size(M,1)
    fprintf('%d', M(i,1));
    
    for j = 2:size(M,2)
        fprintf(' %d', M(i,j));
    end
    
    fprintf('\n');
end