% delete once done!! 205-06-20
function [uop_bps, uop_logprob_persec, crm_bps, crm_logprob_persec, opt_params] = rastbps_comp_findPS(inp_rast,t_bin,ps_basis,options)
% AKHEITMAN: STAND ALONE
% Compute conditioned model (find Post spikefilter) and optimal rate model from inp_rast
% Compute corresponding Bits Per Spike
%{
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
% AKHEITMAN 2015-03-02

%}

% AKHEITMAN 2015-03-01
% Version 0   Works in correcting abeerant WN vs NSEM comparison 
% Version 1_0 Make Everything Stand Alone, import ps_basis
% Version 1_1 Make a Verbose Option
% Version 1_2 Match outputs to rastbps_comp_importPS
% Version 1_3 Clean up make sure, not crossvalidated, LOOKS GOOD!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% TEST CODE %%%
%{
% INPUT
clear; clc
load rastbps_comp_findPS_test_1_2.mat
display('test takes about 30 seconds and ~16 iterations')
options.verbose = true;
[uop_bps, uop_logprob_persec, crm_bps, crm_logprob_persec, opt_params] = rastbps_comp_findPS(inp_rast,t_bin,ps_basis,options);
display(sprintf('## uop_bps = %d; crm_bps = %d', uop_bps, crm_bps))
display(sprintf('## uop_bps_noXVAL = %d; crm_bps_noXVAL_Tol5= %d', uop_bps_correct, crm_bps_correct))
figure(1); plot( opt_params.ps_gain); title('PS Gain Change')
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SET COVARIATES
% OPEN UP THE OPTIONS
if exist('options','var')
    if isfield(options,'verbose') && options.verbose;
        verbose = true;
        display('Verbose Version shows convergence iterations')
    end   
end

% PICK PARAMETERS
reps = size(inp_rast,1);
bins = size(inp_rast,2);

% CONSTRUCT THE UOP BY SMOOTHING WITH 10 MSEC GAUSSIANS
sig_val            = 2; % 10 msecs
convolve_dur       = 20; %
convolve_index     = [-convolve_dur*sig_val:convolve_dur*sig_val];
convolve_vec       = normpdf(convolve_index,0,sig_val) / sum(normpdf(convolve_index,0,sig_val) );

% FIND LOG_RATE AND THEN CORRESPONDING COVARIATE [COV_rate,uop];
summed_spikes      = sum(inp_rast,1);
convolved_raw      = conv( (summed_spikes/reps), convolve_vec);
uop                = convolved_raw( (convolve_dur*sig_val+1):(end-convolve_dur*sig_val));
rate_persec  = 1200*uop;
log_rate     = log(rate_persec);
COV_rate     = repmat(log_rate,[1,reps]);

COV_constant = repmat(ones(1,bins),[1, reps]);


% COVARIATES WHICH ARE REP DEPENDENT ( CONDITIONING TERMS )
PS_CONV = cell(reps,1);
for i_rep = 1:reps
    sp_bins         = find(inp_rast(i_rep,:));
    PS_CONV{i_rep}  = loc_convolvespike_basis(sp_bins,ps_basis',bins);
end

total_bins    = reps * bins;
COV_PS        = NaN(size(ps_basis,2),total_bins);
concat_Spikes = NaN(1,total_bins); 
for i_rep = 1:reps
    index_start = (i_rep-1)*bins + 1;
	index_end   =  i_rep * bins;
    COV_PS(:,index_start:index_end) = PS_CONV{i_rep};
    concat_Spikes(index_start:index_end) = inp_rast(i_rep,:);
end
home_spbins        = find(concat_Spikes);
COV         = [COV_constant ; COV_rate ; COV_PS]; 
COV(find(COV<= -36)) = -(36);  %%% exp(-36) is near machine epsilon

%%
ps_length = size(ps_basis,2);

p_init      = [0 1 zeros(1,ps_length)]';
optim_struct = optimset(...
'derivativecheck','off',...
'diagnostics','off',...  % 
'display','off',...  %'iter-detailed',... 
'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
'GradObj','on',...
'largescale','on',...
'Hessian','on',...
'MaxIter',100,... % you may want to change this
'TolFun',10^(-5),...
'TolX',10^(-9));
if exist('verbose','var') && verbose
    optim_struct = optimset(...
    'derivativecheck','off',...
    'diagnostics','off',...  % 
    'display','iter',...  %'iter-detailed',... 
    'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
    'GradObj','on',...
    'largescale','on',...
    'Hessian','on',...
    'MaxIter',100,... % you may want to change this
    'TolFun',10^(-5),...
    'TolX',10^(-9)   );
end

% RUN OPTIMIZATION ON OFFSET, SCALING RATE, PSFILTER
[pstar fstar eflag output] = fminunc(@(p) rastbps_comp_glm_convexopt(p,COV,home_spbins,t_bin),p_init,optim_struct);

%%
lcif            = pstar' * COV;
crm_ratepersec_concat  = exp(lcif);
crm_ratepersec         = NaN(reps,bins);
for i_rep = 1:reps
	index_start = (i_rep-1)*bins + 1;
	index_end   =  i_rep * bins;
	crm_ratepersec(i_rep,:) = crm_ratepersec_concat(index_start:index_end);
end

crm_rateperbin = t_bin*crm_ratepersec;
crm_model      = crm_rateperbin;
uop_model      = repmat(uop,[reps,1]);
null_model     = repmat( (1/(reps*bins)*sum(inp_rast(:)))*ones(1,bins) , [reps,1]);

% Evaluate Logarthmic Probability
crm_logprob          = sum(rastbps_comp_logprob( inp_rast, crm_model,  'binary')) ;
uop_logprob          = sum(rastbps_comp_logprob( inp_rast, uop_model,  'binary')) ;
null_logprob         = sum(rastbps_comp_logprob( inp_rast, null_model,  'binary')) ;

% Bits Per Spike (improvement from null model)
uop_bps  = (uop_logprob - null_logprob) / (sum(null_model(1,:)));
crm_bps  = (crm_logprob - null_logprob) / (sum(null_model(1,:)));

% LogProb_Sec (raw monotonic measurement of probability of spike train)
uop_logprob_persec = uop_logprob / (bins*t_bin);
crm_logprob_persec = crm_logprob / (bins*t_bin);



opt_params.offset      = pstar(1);
opt_params.rate_drive  = pstar(2);
opt_params.ps_filter   = ps_basis *pstar( (end-ps_length+1):end);
opt_params.ps_gain     = exp(ps_basis *pstar( (end-ps_length+1):end));
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



function [spikesconvbasis]  = loc_convolvespike_basis(binned_spikes,basis,bins)
% AKHeitman 2014-05-04
% Parameter independent!
% basis should be a vector of [basis vectors , bins] already binned
% t_bin is used to put spike times into their appropriate bin 
% t_bin is duration of each bin in msecs
% bins  is the maximum bin number


vectors = size(basis,1); vectorbins= size(basis,2);
%offset by 1,  so PS filter starts in the next bin!!! 
binned_spikes = binned_spikes(find(binned_spikes < bins) ) + 1;


convolvedspikes_base                  = zeros(1,bins);
convolvedspikes_base(binned_spikes+1) = 1; 

convolvedspikes = zeros(vectors, bins + vectorbins - 1);
for i_vec = 1:vectors
    convolvedspikes(i_vec, :) = conv(convolvedspikes_base, basis(i_vec,:), 'full');
end
convolvedspikes = convolvedspikes(:, 1:bins);    

spikesconvbasis = convolvedspikes;
end






    

