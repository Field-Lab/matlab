function [Hess_Correction] = glm_nonconvex_HessianCorrection(log_cif,spikebins,X_frame,frame_shifts,bins_per_frame, bin_duration)
%%% PURPOSE %%%
% Correction term in Hessian due to non-convex parameter use
% Do this in Frames

%%% INPUTS  %%%
% Log_Cif
% Spikebins
% X_frame
% frame_shifts
% bpf

%%% OUTPUTS %%%
% Hess_Correction

% AKHEITMAN 2014-12-15  Version 0

% Version 1 begins in 2015

cif = exp(log_cif);
frames = size(X_frame,2);
cif_byframe_size = reshape(cif,[bins_per_frame,frames]);
cif_byframe      = mean(cif_byframe_size,1);


spikeframes = ceil(spikebins/(bins_per_frame));
spiketime_vec = zeros(1,frames);
for i_spike = 1:length(spikeframes)
    index = spikeframes(i_spike);
    spiketime_vec(index) = spiketime_vec(index) + 1;
end

multiplier_vec_base = (spiketime_vec - bin_duration * bins_per_frame * cif_byframe)';



% Shift the multiplier vec rather than go backwards in time
Hess0 = zeros(size(X_frame,1) , length(frame_shifts));
for i_lag = 1:length(frame_shifts)
    shift = frame_shifts(i_lag);
    multiplier_vec = circshift(multiplier_vec_base, shift);
    Hess0(:,i_lag) = X_frame * multiplier_vec;
end

    
% All computations have this extra multiplier due to fminunc
Hess_Correction = -Hess0;


end