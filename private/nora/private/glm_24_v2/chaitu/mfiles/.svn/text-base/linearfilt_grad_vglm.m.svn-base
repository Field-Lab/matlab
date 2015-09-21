% This function computes the gradient of the input term at each time wrt
% the spatial/temporal filters. The index of the neuron in the param vector
% is specified by offset.

% It returns a (pars.n*pars.Mk) x pars.maxt matrix

% IMPORTANT:
% All of this should be done at the coarser stimulus temporal resolution
% (use stimpars.dt!)

% This version is for nonseparable linear filter

% based on linearfilt_grad_nonsep

function L = linearfilt_grad_vglm(basepars,stimpars)

[n,T] = size(stimpars.x);
% note: for now, i'm assuming that n=1 (edoi)
L = basepars.padval.*ones(n,T); % t-th column is the gradient of the input term at time t
% basepars.padval should be set for vGLM (for the stimulus, it's 0.5)

%-- for the basepars.Mk frame and the rest
%L(basepars.Mk:end) = stim_design(stimpars.x,1);
L(:) = stim_design(stimpars.x,1);

%-- for the first basepars.Mk-1 frames
% for j=1:basepars.Mk-1 
%    %L(1:basepars.n*j,j) = reshape(stimpars.x(space_idx,j:-1:1),j*basepars.n,1);
%    L(j) = reshape(stimpars.x(j:-1:1), j, 1);
% end
