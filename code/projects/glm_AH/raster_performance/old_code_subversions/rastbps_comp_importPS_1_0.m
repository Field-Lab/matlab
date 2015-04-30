function [uop_bps crm_bps] = rastbps_comp(raster,t_bin,psfilter,options)
% STAND-ALONE CODE
% Compute conditioned model and optimal rate model from raster
% Compute corresponding Bits Per Spike

%%% PURPOSE %%%
% Compute a Conditioned Rate Model (CRM) from the raster and psfilter
% Compute the Uncoditioned Optimal Rate Model from raster alone
% Compute corresponding Bits_Per_Spike Code
%
%%% NOTES %%%
% Implement 10 msec guassian smoothing time for all spikes
% Uses post-spike filter generated previously (glm fits) as means for
% conditioning
% Fits run through GLM like fit procedure
% Base Model features only 3 free parameters
%
%%% INPUTS  %%%
% raster: rows are repitions, columns times, binary 0 1 for spike 
% t_bin: time duration of each bin in seconds
% psfilter: log(gain change) induced by a spike
%           concept directly from the GLM
%
%%% OUTPUTS %%%
% uop_bps: Bits per Spike of the Unconditioned Optimal Model
% crm_bps: Bits per Spike of the Conditioned Rate Model
%
%%% OUTSIDE CALLS %%%
% NONE
%
%%% KEY MATLAB CALLS %%%
% fminunc
%
% AKHEITMAN 2015-03-01
% Version 0_0  confirmed to work 2015-03-01
% Version 0_1  Full standalone code: confirmated to work 2015-03-01
% Version 1    Includes Test Code for base case (no xval, simple crm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% TEST CODE %%%
%{
% INPUT
clear; clc
load rastbps_comp_test.mat
[uop_bps crm_bps] = rastbps_comp(raster,t_bin,psfilter);
display('#### TESTING BASE CASE ####')
display(sprintf('## uop_bps = %d; crm_bps = %d', uop_bps, crm_bps))
display(sprintf('## uop_bps_correct = %d; crm_bps_correct = %d', uop_bps_correct, crm_bps_correct))
if uop_bps == uop_bps_correct && crm_bps == crm_bps_correct
    display('## FOR BASE CASE: Exact Same Answer---things look good!')
else
    display('## FOR BASE CASE: Answers not same---check for error!')
end



% Next Test
% DESIRED OUTPUT
% uop_bps = 2.895830e-01; crm_bps = 4.697288e-01
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% PICK PARAMETERS
reps = size(raster,1);
bins = size(raster,2);


% CONSTRUCT THE UOP BY SMOOTHING WITH 10 MSEC GAUSSIANS
binned_rate        = sum(raster,1);
sig_val            = 12; % 10 msecs
convolve_index     = [-1000*sig_val:1000*sig_val];
convolve_vec       = normpdf(convolve_index,0,sig_val) / sum(normpdf(convolve_index,0,sig_val) );
prob_fulluop_bin   = sum(raster,1) / (size(raster,1));
convolved_raw      = conv(prob_fulluop_bin, convolve_vec);
fulluop_smooth     = convolved_raw( (1000*sig_val+1):(end-1000*sig_val));
uop                = fulluop_smooth;  
clear fulluop_smooth convovled_raw prob_fulluop_bin 
clear convolve_vec convolve_index sig_val binned_rate

% CONSTRUCT PS_CONV
PS_CONV = NaN(size(raster));
for i_rep = 1:reps
    sptrain = raster(i_rep,:);
    ps_conv_raw = conv(sptrain,[0;psfilter]);
	PS_CONV(i_rep,:) = ps_conv_raw(1:bins);
end
clear sptrain ps_conv_raw

rate_persec  = 1200*uop;
meanrate     = mean(rate_persec);
gainchange   = rate_persec / meanrate;
log_gain     = log(gainchange);
log_mean     = log(meanrate);
%clear rate_persec meanrate gainchange

% CONSTRUCT THE COVARIATE VECTORS CONSTRUCT
total_bins = reps * bins;
COV_constant  = NaN(1,total_bins);
COV_rate      = NaN(1,total_bins);
COV_postspike = NaN(1,total_bins);
concat_Spikes = NaN(1,total_bins); 
for i_rep = 1:reps
	index_start = (i_rep-1)*bins + 1;
	index_end   =  i_rep * bins;
    concat_Spikes(index_start:index_end) = raster(i_rep,:);
    COV_postspike(index_start:index_end) = PS_CONV(i_rep,:);
    COV_constant(index_start:index_end)  = ones(1,bins);
    COV_rate(index_start:index_end)      = log_gain;    
end

home_spbins = find(concat_Spikes);
COV         = [COV_constant ; COV_rate ; COV_postspike];
p_init = [log_mean 1 0]';

optim_struct = optimset(...
'derivativecheck','off',...
'diagnostics','off',...  % 
'display','off',...  %'iter-detailed',... 
'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
'GradObj','on',...
'largescale','on',...
'Hessian','on',...
'MaxIter',100,... % you may want to change this
'TolFun',10^(-(5)),...
'TolX',10^(-(9))   );

[pstar fstar eflag output] = fminunc(@(p) rastbps_comp_glm_convexopt(p,COV,home_spbins,t_bin),p_init,optim_struct);
lcif            = pstar' * COV;
crm_ratepersec_concat  = exp(lcif);
crm_ratepersec = NaN(reps,bins);
for i_rep = 1:reps
	index_start = (i_rep-1)*bins + 1;
	index_end   =  i_rep * bins;
	crm_ratepersec(i_rep,:) = crm_ratepersec_concat(index_start:index_end);
end

crm_rateperbin = t_bin*crm_ratepersec;
crm_model      = crm_rateperbin;
uop_model      = repmat(uop,[reps,1]);
null_model     = repmat( (1/(reps*bins)*sum(raster(:)))*ones(1,bins) , [reps,1]);

crm_logprob          = sum(rastbps_comp_logprob( raster, crm_model,  'binary')) ;
uop_logprob          = sum(rastbps_comp_logprob( raster, uop_model,  'binary')) ;
null_logprob         = sum(rastbps_comp_logprob( raster, null_model,  'binary')) ;

uop_bits = uop_logprob - null_logprob;
uop_bps  = uop_bits / (sum(null_model(1,:)));

crm_bits = crm_logprob - null_logprob;
crm_bps  = crm_bits / (sum(null_model(1,:)));

end


function  [f grad Hess log_cif]  = rastbps_comp_glm_convexopt(linear_params,covariates,spikebins,bin_duration)
% Brought local copy of optimization function to keep the code stand alone
%
%%%%%%%%% PURPOSE %%%%%%%
% Compute the Objective Function being optimized (f)
% Compute the grad/hess as well for the optimizer
% Monotonically related to  negative of log_cif
% log_cif:= log of the conditional intensity function


%%% NOTES %%%%
% Row indexes the parameters and corresponding covariate
% Column indexes time

%%% INPUTS  %%%
% Params: glm parameters to be optimized, column vector
% Covariates: time dependent input with multiplies into the params
% SpikeBins: spike time in bins of cell
% Bin_Duration: duration of each bin in seconds 
% AKHEITMAN 2014-12-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize
p = linear_params;
COV = covariates;
dt = bin_duration;
spt = spikebins;


% Find Conditional Intensity and its log
lcif = p' * COV;
cif  = exp(lcif);


% Evaluate the objective function (monotonic in log-likelihood)
f_eval = sum( lcif(spt) ) - dt * sum(cif);

% Evaluate the gradient
g_eval = sum(COV(:,spt),2)  - dt * ( COV * (cif') );

% Evaluate the hessian
hessbase = zeros(size(COV));
for i_vec = 1:size(COV,1)
    hessbase(i_vec,:) = sqrt(cif) .* COV(i_vec,:) ;
end
H_eval = -dt * (hessbase * hessbase');


% Switch signs because using a minimizer  fmin
f       = -f_eval;
grad    = -g_eval;
Hess    = -H_eval;
log_cif = lcif;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ raster_logprob_bin, logprobmat] = rastbps_comp_logprob( spikemat, ratemodel, spikemat_type) 
% Brought local copy of optimization function to keep the code stand alone
% 
%%% PURPOSE %%%
% Compute Likelihood of Spike Raster (usa binary)
% Model is a binned rate model

%%% NOTES %%%%
% Row Indexes the repitition
% Column indexes time


trials = size(spikemat,1);
bins   = size(spikemat,2);

if strcmp(spikemat_type, 'binary')
    % EVALUATE MODEL LIKELIHOOD OF SPIKE AND NO SPIKES
    P_nospike = exp(-ratemodel);
    P_spike   = ratemodel.*P_nospike;
    
    % FIND SPIKES
    spindex   = find(spikemat);
    
    % EVALUATE PROBABLITY OF SPIKEMAT GIVEN MODEL
    probmat             = P_nospike;
    probmat(spindex)    = P_spike(spindex); 
    logprobmat          = log(probmat);
    raster_logprob_bin  = mean(logprobmat,1);
end

end



