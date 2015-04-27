% Compute the convolution of the linear temporal filter basis fns with the
% stimulus. This is used in computing the hessian wrt space/time filters.

function [XP] = hess_spacetimeconst(x,phi,dt)

[n t] = size(x);

[Mk N] = size(phi);

XP = zeros(n,t,N); % the kth slice is stimulus convolved with kth temporal filter

for k=1:N
    % Convolve the kth basis fn with the stimulus - take only the first t cols
    XP(:,:,k) = dt*fastconv(x,phi(:,k)',n,t);
end