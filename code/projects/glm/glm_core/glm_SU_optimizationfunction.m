function  [f grad Hess log_cif] = glm_SU_optimizationfunction(SU_params, SU_covariates, pooling_weights, spikebins, bin_duration, non_stim_lcif)
 %%% PURPOSE %%%
% Compute the Objective Function being optimized (f)
% Compute the grad/hess as well for the optimizer
% Monotonically related to  negative of log_cif
% log_cif:= log of the conditional intensity function


%%% This one is designed for each linear_param to be a subunit weight (eg
%%% go through a nonlinearity)

%%% NOTES %%%%
% Row indexes the parameters and corresponding covariate
% Column indexes time

%%% INPUTS  %%%
% Params: glm parameters to be optimized, column vector
% Covariates: time dependent input with multiplies into the params
% SpikeBins: spike time in bins of cell
% Bin_Duration: duration of each bin in seconds 


% AKHEITMAN 2014-12-04
%%

% Initialize
dt = bin_duration;
spt = spikebins;
SU_covariates = imresize(SU_covariates, [1, length(non_stim_lcif)],'nearest');

% hessbase = zeros(size(COV));

% Find Conditional Intensity and its log
pixels = SU_covariates;
for i = 1:length(SU_params)
    SU_covariates(i, :, :) = SU_covariates(i, :, :) * SU_params(i);
end
stim_lcif = pooling_weights'* exp(squeeze(sum(SU_covariates,1)));
lcif = stim_lcif + non_stim_lcif;

% Evaluate the objective function (monotonic in log-likelihood)
f_eval = sum( lcif(spt) ) - dt * sum(exp(lcif));
% 
% % Evaluate the gradient

% gradient of the lcif
del_lcif = zeros(9, length(SU_covariates));
subunit_drive = exp(squeeze(sum(SU_covariates,1)));
for i_SU = 1:9
    for i_loc = 1:121
        temp = pooling_weights(i_loc)*pixels(i_SU,i_loc,:).*subunit_drive(i_loc,:);
        del_lcif(i_SU,:) = del_cif(i_SU,:)+temp;
    end
end
g_eval = sum(COV(:,spt),2)  - dt * ( COV * (cif') );
% 
% % Evaluate the hessian
% hessbase(i_vec,:) = sqrt(cif) .* COV(i_vec,:) ;
% 
% H_eval = -dt * (hessbase * hessbase');
% 


%%

% Switch signs because using a minimizer  fmin
f       = -f_eval;
% grad    = -g_eval;
% Hess    = -H_eval;
% log_cif = lcif;
end