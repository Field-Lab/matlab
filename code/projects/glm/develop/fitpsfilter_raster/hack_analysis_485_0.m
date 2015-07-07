% 2015-06-22
% Begin exploration of Post-Spike filter
% Simulate.. pull raster.. see what happens.
% Really a test of concept of rastbps_wrap


clear; close all
load('cell485_fitfromrast_NSEM.mat')
load('cell485_baserate_fromrast.mat')

sim.log_input_rate = opt_params.rate_drive*log(rate_persec)+opt_params.offset
sim.ps_gain        = opt_params.ps_gain;

params.trials = 2*size(inp_rast,1);
params.bins   = size(inp_rast,2);
params.bindur = 1/1200;   % hack in
t_bin = 1/1200;



% Unpack inputs
% lcif_rate = repmat(sim.log_input_rate , [params.trials, 1]); 
lcif_rate  = sim.log_input_rate;
cif_psgain = sim.ps_gain;
ps_bins    = length(cif_psgain);
% Simulate
logical_sim = zeros(params.trials, params.bins);
for i_trial = 1:params.trials
    display(sprintf('simulating trial %d',i_trial))
    cif0         = exp(lcif_rate);
    cif_ps       = cif0;
    binary_simulation = zeros(1,params.bins);
    for i = 1 : params.bins- ps_bins;
        roll = rand(1);
        if roll >  exp(-params.bindur*cif_ps(i));
            cif_ps(i+1: i + ps_bins) =  cif_ps(i+1: i + ps_bins) .* (cif_psgain');
            binary_simulation(i)= 1;
        end
    end
    logical_sim(i_trial,:) = binary_simulation ;
end


rast = logical_sim;
% CONSTRUCT THE UOP BY SMOOTHING WITH 10 MSEC GAUSSIANS
reps = size(rast,1);
sig_val            = 12; % 10 msecs
convolve_dur       = 10; %
convolve_index     = [-convolve_dur*sig_val:convolve_dur*sig_val];
convolve_vec       = normpdf(convolve_index,0,sig_val) / sum(normpdf(convolve_index,0,sig_val) );

% FIND LOG_RATE AND THEN CORRESPONDING COVARIATE [COV_rate,uop];
summed_spikes      = sum(rast,1);
convolved_raw      = conv( (summed_spikes/reps), convolve_vec);
uop                = convolved_raw( (convolve_dur*sig_val+1):(end-convolve_dur*sig_val));
newrate_persec  = 1200*uop;



tic
[uop_bps, uop_logprob_persec, crm_bps, crm_logprob_persec, opt_params2] = rastbps_comp_findPS(logical_sim,t_bin,ps_basis);
toc

sim2.log_input_rate = opt_params2.rate_drive*log(newrate_persec)+opt_params2.offset;
sim2.ps_gain  = opt_params2.ps_gain;









% Unpack inputs
% lcif_rate = repmat(sim2.log_input_rate , [params.trials, 1]); 
lcif_rate  = sim2.log_input_rate;
cif_psgain = sim2.ps_gain;
ps_bins    = length(cif_psgain);
% Simulate
logical_sim = zeros(params.trials, params.bins);
for i_trial = 1:params.trials
    display(sprintf('simulating trial %d',i_trial))
    cif0         = exp(lcif_rate);
    cif_ps       = cif0;
    binary_simulation = zeros(1,params.bins);
    for i = 1 : params.bins- ps_bins;
        roll = rand(1);
        if roll >  exp(-params.bindur*cif_ps(i));
            cif_ps(i+1: i + ps_bins) =  cif_ps(i+1: i + ps_bins) .* (cif_psgain');
            binary_simulation(i)= 1;
        end
    end
    logical_sim(i_trial,:) = binary_simulation ;
end


rast = logical_sim;
% CONSTRUCT THE UOP BY SMOOTHING WITH 10 MSEC GAUSSIANS
reps = size(rast,1);
sig_val            = 12; % 10 msecs
convolve_dur       = 10; %
convolve_index     = [-convolve_dur*sig_val:convolve_dur*sig_val];
convolve_vec       = normpdf(convolve_index,0,sig_val) / sum(normpdf(convolve_index,0,sig_val) );

% FIND LOG_RATE AND THEN CORRESPONDING COVARIATE [COV_rate,uop];
summed_spikes      = sum(rast,1);
convolved_raw      = conv( (summed_spikes/reps), convolve_vec);
uop                = convolved_raw( (convolve_dur*sig_val+1):(end-convolve_dur*sig_val));
newnewrate_persec  = 1200*uop;

clf
plot(rate_persec,'b'); hold on; plot(newrate_persec,'r');  plot(newnewrate_persec,'m'); 
%plot(exp(sim.log_input_rate), 'm')

[uop_bps, uop_logprob_persec, crm_bps, crm_logprob_persec, opt_params3] = rastbps_comp_findPS(logical_sim,t_bin,ps_basis);
figure; 
plot(opt_params.ps_gain,'b'); hold on;plot(opt_params2.ps_gain,'r'); plot(opt_params3.ps_gain,'m'); 

[uop_bps, uop_logprob_persec, crm_bps, crm_logprob_persec, opt_params_shortsmooth] = rastbps_comp_findPS_shortsmooth(logical_sim,t_bin,ps_basis);
[uop_bps, uop_logprob_persec, crm_bps, crm_logprob_persec, opt_params_longsmooth] = rastbps_comp_findPS_longsmooth(logical_sim,t_bin,ps_basis);
figure;
hold on;plot(opt_params_shortsmooth.ps_gain,'k'); plot(opt_params3.ps_gain,'m');  plot(opt_params_longsmooth.ps_gain,'r'); 