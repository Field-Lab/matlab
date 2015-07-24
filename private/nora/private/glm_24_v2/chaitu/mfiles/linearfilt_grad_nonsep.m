% This function computes the gradient of the input term at each time wrt
% the spatial/temporal filters. The index of the neuron in the param vector
% is specified by offset.

% It returns a (pars.n*pars.Mk) x pars.maxt matrix

% IMPORTANT:
% All of this should be done at the coarser stimulus temporal resolution
% (use stimpars.dt!)

% This version is for nonseparable linear filter

%function L = linearfilt_grad_nonsep(p,basepars,stimpars,trainpars,neuron_no)
function L = linearfilt_grad_nonsep(~,basepars,stimpars,~,neuron_no)

[~, T] = size(stimpars.x);
nspace = basepars.n;

if (iscell(basepars.crop_idx))
   space_idx = basepars.crop_idx{neuron_no};
else
   space_idx = basepars.crop_idx(:,neuron_no);
end

L = basepars.padval.*ones(nspace*basepars.Mk,T); % t-th column is the gradient of the input term at time t

%-- for the first basepars.Mk-1 frames
for j=1:basepars.Mk-1 
   L(1:basepars.n*j,j) = reshape(stimpars.x(space_idx,j:-1:1),j*basepars.n,1);
end
%-- for the basepars.Mk frame and the rest
L(:,basepars.Mk:end) = stim_design(stimpars.x(space_idx,:),basepars.Mk);

%-- not used (edoi)
if (isfield(basepars,'rms_mlength') && basepars.rms_mlength > 0)
   rms = compute_rms(L,basepars.rms_mlength,basepars.rms_slength,basepars.padval);
   L = L./rms;
end
