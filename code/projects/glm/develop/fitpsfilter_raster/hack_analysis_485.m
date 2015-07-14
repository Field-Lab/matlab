% 2015-06-22
% Begin exploration of Post-Spike filter
% Simulate.. pull raster.. see what happens.
% Really a test of concept of rastbps_wrap


clear; close all
load('cell485_fitfromrast_NSEM.mat')
load('cell485_baserate_fromrast.mat')

sim.log_input_rate = opt_params.rate_drive*log(rate_persec)+opt_params.offset
sim.ps_gain        = opt_params.ps_gain;

opt_params_sim = opt_params; clear opt_params;
%% showing not a problem of smoothing
options.verbose = true;

smooth_rates = [1 3 6 12 24 48 96];
results = cell(length(smooth_rates),1);
for i_rate = 1:length(smooth_rates)
    % should roughly take a minute each
    [uop_bps, ~,crm_bps, ~, opt_params] = rastbps_comp_findPS_variabletime(inp_rast,...
        1/1200,ps_basis,smooth_rates(i_rate),options);
    results{i_rate}.uop_bps = uop_bps;
    results{i_rate}.crm_bps = crm_bps;
    results{i_rate}.opt_params = opt_params;
    
    clear opt_params;
end

ps_basis_unit


figure;
set(gca, 'fontsize',12); hold on
smooth_rates = [1 3 6 12 24 48 96];
lw_other = 1.2;
time = linspace(0,100,120);
for i_rate = 1:length(smooth_rates)
    if i_rate == 1; color = 'r'; end
    if i_rate == 2; color = 'r'; end
    if i_rate == 3; color = 'r'; end
    if i_rate == 4; color = 'k'; end
    if i_rate == 5; color = 'b'; end
    if i_rate == 6; color = 'g'; end
    if i_rate == 7; color = 'm'; end
    
    if i_rate == 4
        plot(time,results{i_rate}.opt_params.ps_gain,color,'linewidth', 3); hold on
    else
        plot(time,results{i_rate}.opt_params.ps_gain,color,'linewidth', lw_other); hold on
    end
end
xlabel('msecs'); ylabel('gain');
xlim([0 100]); ylim([0 1.25]);
legend('1 msec','2.5msec','5msec','10msec(defaultscale)','20msec','40msec','80msec','location','southeast')
title('PS filter fitted from 10msec smoothed raster is not too smooth')
plot(time,ones(size(time)), 'k')      

figure
set(gca, 'fontsize',12); hold on
smooth_rates = [1 3 6 12 24 48 96];
lw_other = 1.1;
time = linspace(0,100,120);
for i_rate = 1:length(smooth_rates)
    if i_rate == 1; color = 'r'; end
    if i_rate == 2; color = 'r'; end
    if i_rate == 3; color = 'r'; end
    if i_rate == 4; color = 'k'; end
    if i_rate == 5; color = 'b'; end
    if i_rate == 6; color = 'g'; end
    if i_rate == 7; color = 'm'; end
    
    if i_rate == 4
        plot(time,results{i_rate}.opt_params.ps_gain,color,'linewidth', 3); hold on
    else
        plot(time,results{i_rate}.opt_params.ps_gain,color,'linewidth', lw_other); hold on
    end
end
xlabel('msecs'); ylabel('gain');
xlim([0 50]); ylim([0 1.25]);
legend('1 msec','2.5msec','5msec','10msec(defaultscale)','20msec','40msec','80msec','location','southeast')
title('PS filter fitted from 10msec Smoothed Raster ')
plot(time,ones(size(time)), 'k')   

if strcmp(crm_type, 'base_crm_findPS') || strcmp(crm_type, 'base_crm_findPS_HO')
    % HARD PS PARAMETERS  MATCHES GLMPARAMS
    t_bin0 = .00083275;
    ps_params.ms           = 100 ;     
    ps_params.filternumber = 20;
    ps_params.spacing      = pi/2;
    ps_params.bstretch     = .05;
    ps_params.alpha        = 0;
    ps_params.fratio       = .5;  
    ps_basis               = prep_spikefilterbasisGP(ps_params,t_bin0);
end




%% weaken 
lcif_rate0  = sim.log_input_rate;
lcif_rate1  = log(20+ (exp(lcif_rate0)).^(.5) );
%lcif_rate1  = log(60+ zeros(size(lcif_rate0) ));


sim_subbins = 4;
sim_bin     = 1/1200 * (1/sim_subbins);

params.trials = 2*size(inp_rast,1);
params.bins   = sim_subbins * size(inp_rast,2);
params.bindur = sim_bin;

sim.ps_gain_original = sim.ps_gain;
sim.ps_gain(45:end)  = 1;
sim.ps_gain          = sim.ps_gain(1:60);
ps_basis_restrict    = ps_basis(1:120,1:25);

lcif_rate = repmat(lcif_rate1, [sim_subbins,1]);
lcif_rate = reshape(lcif_rate, 1, params.bins);
cif_psgain= repmat(sim.ps_gain',[sim_subbins,1]);
cif_psgain= reshape(cif_psgain, 1, 60*sim_subbins);
ps_bins    = length(cif_psgain);
% Simulate
logical_sim = zeros(params.trials, params.bins);
for i_trial = 1:params.trials
    if mod(i_trial,10) == 0
        display(sprintf('simulating trial %d out of %d',i_trial,params.trials));
    end
    cif0         = exp(lcif_rate);
    cif_ps       = cif0';
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
newcols  = size(logical_sim,2) / sim_subbins;
new_rast = zeros( size(logical_sim,1) , newcols);
for i_col = 1:newcols
    columns = (i_col-1)*sim_subbins + [1:sim_subbins];
    new_rast(:,i_col) = sum( logical_sim(:,columns),2 );
end
new_rast(find(new_rast>=2)) = 1;
sim_rast = new_rast;
[uop_bps, ~,crm_bps, ~, opt_params_sim] = rastbps_comp_findPS_variabletime(sim_rast,...
        1/1200,ps_basis_restrict,12,options);

figure;

plot(opt_params_sim.ps_gain,'k','linewidth', 1.5)  
hold on; plot(sim.ps_gain,'r','linewidth',1.5);
legend('PS of Sim Input','Fitted PS from Sim Output','location', 'southeast')




%% show process is "robust"  failed attempt, ps filter flattens out too much
%{
sim_subbins = 10;
sim_bin     = 1/1200 * (1/sim_subbins);

params.trials = size(inp_rast,1);
params.bins   = sim_subbins * size(inp_rast,2);
params.bindur = sim_bin;

sim.ps_gain_original = sim.ps_gain;
sim.ps_gain(45:end)  = 1;
sim.ps_gain          = sim.ps_gain(1:60);
ps_basis_restrict    = ps_basis(1:120,1:25);

% Unpack inputs
% lcif_rate = repmat(sim.log_input_rate , [params.trials, 1]); 
lcif_rate0  = sim.log_input_rate;

lcif_rate = repmat(lcif_rate0, [sim_subbins,1]);
lcif_rate = reshape(lcif_rate, 1, params.bins);
cif_psgain= repmat(sim.ps_gain',[sim_subbins,1]);
cif_psgain= reshape(cif_psgain, 1, 60*sim_subbins);
ps_bins    = length(cif_psgain);
% Simulate
logical_sim = zeros(params.trials, params.bins);
for i_trial = 1:params.trials
    if mod(i_trial,10) == 0
        display(sprintf('simulating trial %d out of %d',i_trial,params.trials));
    end
    cif0         = exp(lcif_rate);
    cif_ps       = cif0';
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

newcols  = size(logical_sim,2) / sim_subbins;
new_rast = zeros( size(logical_sim,1) , newcols);
for i_col = 1:newcols
    columns = (i_col-1)*sim_subbins + [1:sim_subbins];
    new_rast(:,i_col) = sum( logical_sim(:,columns),2 );
end
%new_rast(find(new_rast>=2)) = 1;

options.verbose = true;
sim_rast = new_rast;
[uop_bps, ~,crm_bps, ~, opt_params_sim] = rastbps_comp_findPS_variabletime(sim_rast,...
        1/1200,ps_basis_restrict,12,options);

plot(opt_params_sim.ps_gain,'m','linewidth', 1.5)  
hold on; plot(sim.ps_gain,'r')
%}
