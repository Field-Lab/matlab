function [out inds] = isrange(x,mn,mx)

%-------------------------------------------
%
% out = isrange(x,mn,mx)
%
% check whether a set of values are within a given range
%
% input:
% x -- set of values to check (scalar or vector or matrix)
% mn -- lower bound of range
% mx -- upper bound of range
%
% output:
% out -- 1 if all values are within the range, 0 otherwise
% inds -- return the indices
%
% freeman, 5-4-2012
%-------------------------------------------

x = x(:);
inds = x < mx & x > mn;
if sum(inds) == length(x)
    out = logical(1);
else
    out = logical(0);
end