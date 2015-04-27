% Function that does ideal sinc interpolation
% Chaitu Ekanadham

% Arguments:
% taxis - points at which to query function value
% timepts - points at which samples are taken

% Returns:
% A - the interpolation matrix. i.e. if x is a column vector of sample
% values at timepts, then A*x is the (bandlimited) function at taxis,
% assuming zero outside.

function A = sincinterpmtx(taxis,timepts)

if (length(timepts)>1)
    T = timepts(2)-timepts(1);
else
    T = taxis(end);
end

A = sinc((1/T).*(repcols(taxis(:),length(timepts)) - reprows(timepts(:)',length(taxis))));
