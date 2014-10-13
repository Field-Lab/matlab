function [combinations,varyingParamValues] = findCombinations(M,n)
% This function returns an array of cell that contains combinations of
% experiments that should be plotted together when given an experiment
% matrix and a parameter of interest.
% 
% Input:
%   - M: columns of M are parameters, rows of M are experiment run.
%   - n: index of the parameter of interest. 
%
% Output:
%   - an array of cells, that gives the indices of all intesting row
%   combinations for the specified parameter n.
%   - another array of cells, that gives the value of the varying parameter
%   n for the corresponding combination in the first array of cells.
%
% Version: v4.04 - 05/08/2011

subM = M;
subM(:,n) = [];

nRows = size(M,1);
rowsUsed = [];

currentCombination = 1;
combinations = cell(1,1);
varyingParamValues = cell(1,1);

for kk=1:nRows
    if isempty(find(rowsUsed==kk, 1))
        currentRow = subM(kk,:);
        combinations{currentCombination} = [];
        varyingParamValues{currentCombination} = [];
        for ll=kk:nRows
            if  isempty(find(~(subM(ll,:)==currentRow), 1))
                combinations{currentCombination} = [combinations{currentCombination} ll];
                varyingParamValues{currentCombination} = [varyingParamValues{currentCombination} M(ll,n)];
                rowsUsed = [rowsUsed ll];
            end
        end
        currentCombination = currentCombination + 1;
    end
end

end % findCombinations