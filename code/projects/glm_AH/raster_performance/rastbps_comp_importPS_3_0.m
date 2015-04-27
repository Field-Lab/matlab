function [uop_bps uop_logprob_persec crm_bps crm_logprob_persec pstar] = rastbps_comp_importPS(inp_rast,t_bin,psfilter,options)
% AKHEITMAN: STAND-ALONE CODE
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
% Version 1_0  Includes Test Code for base case (no xval, simple crm)
% Version 1_1  Includes Test Code for base case (no xval, simple crm)
% Version 1_2  Hack that gets around 0 firing rate .. 

% Version 2_0 Works, no-tol
% Version 2_1 XVAL and New test.mat .. log-rate rather than log_gain
% Version 2_2 TolX Set to 5, Safegaurd UOP from having 0 rate (machine eps)
% Version 3_0 Added logprob_persec as additional output
%             Make no xval default, xval an option
%             Init point set to [0 1 1] encourage PSfilter adoption
%             TolX set to 6 to guarantee doesn't get stuck
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% TEST CODE %%%
%{
% INPUT
clear; clc
load rastbps_comp_importPS_test_3_0.mat
[uop_bps uop_logprob_persec crm_bps crm_logprob_persec pstar] = rastbps_comp_importPS(inp_rast,t_bin,psfilter);
display('#### TESTING BASE CASE ####')
display(sprintf('## uop_bps = %d; crm_bps = %d', uop_bps, crm_bps))
display(sprintf('## uop_bps_tol6 = %d; crm_bps_tol6 = %d', uop_bps_tol6, crm_bps_tol6))


%clear; clc
%load rastbps_comp_importPS_test_2_0.mat
%[uop_bps crm_bps pstar] = rastbps_comp_importPS(inp_rast,t_bin,psfilter);
%display('#### TESTING BASE CASE ####')
%display(sprintf('## uop_bps = %d; crm_bps = %d', uop_bps, crm_bps))
%display(sprintf('## uop_bps_noxval_tol3 = %d; crm_bps_noxval_tol3 = %d', uop_bps_correct, crm_bps_correct))
%pstar
%if uop_bps == uop_bps_correct && crm_bps == crm_bps_correct 
%    display('## FOR BASE CASE: Exact Same Answer---things look good!')
%else
%    display('## FOR BASE CASE: Answers not same---check for error!')
%end

% Next Test
% DESIRED OUTPUT
% uop_bps = 2.895830e-01; crm_bps = 4.697288e-01
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SET COVARIATES
if exist('options','var')
    if isfield(options, 'XVAL')
        crossval = true; 
    end
end
% PICK PARAMETERS
reps = size(inp_rast,1);
bins = size(inp_rast,2);

% CONSTRUCT THE UOP BY SMOOTHING WITH 10 MSEC GAUSSIANS
sig_val            = 12; % 10 msecs
convolve_dur       = 10; %
convolve_index     = [-convolve_dur*sig_val:convolve_dur*sig_val];
convolve_vec       = normpdf(convolve_index,0,sig_val) / sum(normpdf(convolve_index,0,sig_val) );

% FIND LOG_RATE AND THEN CORRESPONDING COVARIATE�[COV_rate,uop];
summed_spikes      = sum(inp_rast,1);
convolved_raw      = conv( (summed_spikes/reps), convolve_vec);
uop                = convolved_raw( (convolve_dur*sig_val+1):(end-convolve_dur*sig_val));
rate_persec  = 1200*uop;
log_rate     = log(rate_persec);
COV_rate     = repmat(log_rate,[1,reps]);
clear log_rate rate_persec convolved_raw summed_spikes
% Optional: exercise "crossvalidated" rate
if exist('crossval','var') && crossval
    summed_spikes      = sum(inp_rast,1);
    prob_spbin_XVAL    = NaN(1,reps*bins);
    for i_rep = 1:reps
        index_start = (i_rep-1)*bins + 1;
        index_end   =  i_rep * bins;
        new_term    = (summed_spikes - inp_rast(i_rep,:)) / (reps -1 ) ;
        prob_spbin_XVAL(index_start:index_end) = new_term;
    end
    convolved_raw      = conv(prob_spbin_XVAL, convolve_vec);
    fulluop_smooth     = convolved_raw( (convolve_dur*sig_val+1):(end-convolve_dur*sig_val));
    uop                = fulluop_smooth;  
    rate_persec  = 1200*uop;
    log_rate     = log(rate_persec);
    COV_rate     = log_rate;
    clear log_rate rate_persec uop convolved_raw summed_spikes
    clear fulluop_smooth prob_sp_bin_XVAL index_start index_end
end


% CONSTRUCT PS_CONV
PS_CONV = NaN(size(inp_rast));
for i_rep = 1:reps
    sptrain = inp_rast(i_rep,:);
    ps_conv_raw = conv(sptrain,[0;psfilter]);
	PS_CONV(i_rep,:) = ps_conv_raw(1:bins);
end


% COVARIATES VWHICH ARE REP INDEPENDENT (RATE, MEAN)
COV_constant   = repmat(ones(1,bins),[1, reps]);

% COVARIATES WHICH ARE REP DEPENDENT ( CONDITIONING TERMS )
total_bins    = reps * bins;
COV_postspike = NaN(1,total_bins);
concat_Spikes = NaN(1,total_bins); 
for i_rep = 1:reps
	index_start = (i_rep-1)*bins + 1;
	index_end   =  i_rep * bins;
    concat_Spikes(index_start:index_end) = inp_rast(i_rep,:);
    COV_postspike(index_start:index_end) = PS_CONV(i_rep,:); 
end

home_spbins = find(concat_Spikes);
COV         = [COV_constant ; COV_rate ; COV_postspike] ;

if exist('higher_powers','var') && higher_powers
    COV         = [COV_constant ; COV_rate ; COV_postspike; COV_rate_p2; COV_rate_n2; COV_rate_p3; COV_rate_n3];
end
% HACK!  Set lowest log_gain to log(machine_eps) 
% Note   log(machine_eps) = -36
COV(find(COV<= -36)) = -(36);

%% OPTIMIZE THE CONDITIONED MODEL
% INITIALIZE AND OPTIMIZE
p_init  = [0 1 1]';
optim_struct = optimset(...
'derivativecheck','off',...
'diagnostics','off',...  % 
'display','off',...  %'iter-detailed',... 
'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
'GradObj','on',...
'largescale','on',...
'Hessian','on',...
'MaxIter',100,... % you may want to change this
'TolFun',10^(-(6)),...
'TolX',10^(-(9))   );
[pstar fstar eflag output] = fminunc(@(p) rastbps_comp_glm_convexopt(p,COV,home_spbins,t_bin),p_init,optim_struct);


%% EVALUATE BOTH CONDITIONED MODEL AND THE PURE RATE MODEL
% DEFINE THE CRM_MODEL
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

% DEFINE THE UOP_MODEL
uop_model      = repmat(uop,[reps,1]);
if exist('crossval','var') && crossval
    uop_model_concat = uop;
    uop_model        = NaN(reps,bins);
    for i_rep = 1:reps
        index_start = (i_rep-1)*bins + 1;
        index_end   =  i_rep * bins;
        uop_model(i_rep,:) = uop_model_concat(index_start:index_end);
    end
    uop_model(find(uop_model<=0)) = eps; 
end

% Define Null Model (optimal tonic firing rate)
null_model           = repmat( (1/(reps*bins)*sum(inp_rast(:)))*ones(1,bins) , [reps,1]);

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
end
%% CALLED SUBFUNCTIONS

function  [f grad Hess log_cif] = rastbps_comp_glm_convexopt(linear_params,covariates,spikebins,bin_duration)
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



