% AKHeitman 2015-07-02
% Just refit PS

%{
clear; close all
%load('dbug_glmexecute_constrainPS_debugWNfit.mat')
load('dbug_glmexecute_constrainPS_fullWNfit.mat')
domainconstrain_name = 'PS_inhibitorydomainconstrain_post10msec'; 
[fittedGLM] = glm_execute_domainconstrainPS(domainconstrain_name, GLMType,...
    fitspikes_concat,fitmovie_concat,testspikes_raster,testmovie,...
    inputstats,glm_cellinfo); 
%}
function [fittedGLM] = glm_execute_refitPS(newPS, GLMType,stimdrivenrate,fitspikes,testspikes_logicalraster,glm_cellinfo)


%% Setup Covariates
fittedGLM.cell_savename = glm_cellinfo.cell_savename;
fittedGLM.d_save        = glm_cellinfo.d_save;
fittedGLM.cellinfo      = glm_cellinfo;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load up GLMParams compute some universal params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GLMPars           = GLMParams;
if isfield(GLMType, 'specialchange') && GLMType.specialchange
    GLMPars = GLMParams(GLMType.specialchange_name);
end
% IF RUNNING OFF OF INPUT NL LOAD THOSE COEFFS
if isfield(glm_cellinfo, 'GLMPars')
   GLMPars = glm_cellinfo.GLMPars; 
end

% Perhaps we should combine this! With convolving with spikes !
t_bin      = glm_cellinfo.t_bin;
if GLMType.PostSpikeFilter
    basis_params  = GLMPars.spikefilters.ps;
    ps_basis      = prep_spikefilterbasisGP(basis_params,t_bin);
end
if GLMType.CouplingFilters
    basis_params  = GLMPars.spikefilters.cp;
    cp_basis      = prep_spikefilterbasisGP(basis_params,t_bin);
end

if strcmp(newPS,'netinhibitory_domainconstrain_COB')
    ps_basis_0 = ps_basis; clear ps_basis
    v        = sum(ps_basis_0,1);
    v        = v / norm(v) ;
    orthog_v = null(v);
    COB      = [v', orthog_v] ;
    ps_basis = (inv(COB) * ps_basis_0')' ;
    clear v orthog_v COB
end

fittedGLM.GLMPars = GLMPars;
fittedGLM.GLMType = GLMType;
if isfield(GLMType, 'debug') && GLMType.debug
    GLMPars.optimization.tolfun = 1; 
end

% Convolve Spike Times with appropriate basis
% Think about flushing dt out to the wrapper
% Take care of all timing in glm_execute or in glmwrap.
bins = length(stimdrivenrate.fit)
home_sptimes = fitspikes.home';
home_spbins  = ceil(home_sptimes / t_bin);
home_spbins = home_spbins(find(home_spbins < bins) );


basis         = ps_basis';
PS_bin        = prep_convolvespikes_basis(home_spbins,basis,bins);
MU_bin   = ones(1,bins);
STIM_bin = log(stimdrivenrate.fit);

glm_covariate_vec = [MU_bin; STIM_bin; PS_bin];

paramcount = size(glm_covariate_vec,1);

%%
% PREPARE PARAMETERS

if strcmp(newPS,'netinhibitory_domainconstrain_COB') 
    lowerbound = -Inf(paramcount,1);
    upperbound  = Inf(paramcount,1);
    upperbound(3) = 0;
    
    fittedGLM.constrained_search.note = 'how the parameter search was limited in fmincon';
    fittedGLM.constrained_search.lowerbound = lowerbound;
    fittedGLM.constrained_search.upperbound = upperbound;
        
    %%%
    optim_struct = optimset(...
   'Algorithm','trust-region-reflective',...
   'derivativecheck','off',...
   'diagnostics','off',...  % 
   'display','iter',...  %'iter-detailed',... 
   'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
   'GradObj','on',...
   'largescale','on',...
   'Hessian','on',...
   'MaxIter',GLMPars.optimization.maxiter,... % you may want to change this
   'TolFun',10^(-(GLMPars.optimization.tolfun)),...
   'TolX',10^(-(GLMPars.optimization.tolx))   ) ;
end


p_init     = [zeros(paramcount,1)];
p_init(2)  = 1;
if isfield(glm_cellinfo, 'p_init')
    p_init = glm_cellinfo.p_init;
    if strcmp(newPS,'netinhibitory_domainconstrain_COB')
        p_init_psbase = p_init(3:end);
        p_init_new    = inv(COB) * p_init_psbase;        
        p_init(3:end) = p_init_new;
        display('modulating p init to new basis')
    end
    
end


%% Run through optimization .. get out pstart, fstar, eflag, output
% CONVEXT OPTIMIZATION
if GLMType.CONVEX
    if strcmp(newPS,'netinhibitory_domainconstrain_COB')
        [pstar fstar eflag output] = fmincon(@(p) ...
            glm_convex_optimizationfunction(p,glm_covariate_vec,home_spbins,t_bin),...
            p_init,[],[],[],[],lowerbound,upperbound,[],optim_struct);
    end
    [f grad Hess log_cif] = glm_convex_optimizationfunction(pstar,glm_covariate_vec,home_spbins,t_bin);
    
    
    
    %{
    [pstar fstar eflag output] = fminunc(@(p) ...
        glm_convex_optimizationfunction(p,glm_covariate_vec,home_spbins,t_bin),...
        p_init,optim_struct);
    [f grad Hess log_cif] = glm_convex_optimizationfunction(pstar,glm_covariate_vec,home_spbins,t_bin);
    plot(ps_basis * pstar_con(paramind.PS),'r'); hold on
    plot(ps_basis * pstar(paramind.PS),'k'); hold on
    % Before constrains
    [pstar fstar eflag output] = fminunc(@(p) ...
        glm_convex_optimizationfunction_constrainPS(p,glm_covariate_vec,home_spbins,t_bin, paramind.PS,ps_basis,PS_balance),...
        p_init,optim_struct);
    [f grad Hess log_cif ps_filter] = glm_convex_optimizationfunction_constrainPS(pstar,glm_covariate_vec,home_spbins,t_bin, paramind.PS,ps_basis,PS_balance);
    %}
end


fittedGLM.fminunc_output = output;

%% Unpack the output into filters

rawfit.opt_params        = pstar;
rawfit.objective_val     = fstar;
rawfit.ps_basis = ps_basis;
% SAVE ALL FILTERS EXCEPT FOR STIMULUS FILTERS
clear linearfilters
linearfilters.note = 'These terms, convolved with the covariates, make the log of the conditional intensity function';
linearfilters.TonicDrive.Filter      = pstar(1);
linearfilters.TonicDrive.note        ='no convolution necessary, this term is part of the lcif for every bin';
linearfilters.stimrescale            = pstar(2);
linearfilters.PostSpike.Filter     = ps_basis * pstar(3:end)
linearfilters.PostSpike.startbin   = 1;  
linearfilters.PostSpike.note0       = 'Filter starts at "startbin" bins after the spikebin';
linearfilters.PostSpike.note0       = 'Filter starts at "startbin" bins after the spikebin';

fittedGLM.rawfit               = rawfit;
fittedGLM.linearfilters = linearfilters;
fittedGLM.note = 'in theory, linearfilters and t_bin/ binsperframe is sufficient for xval and simulation'; 
fittedGLM.fit_time = datestr(clock);
fittedGLM.writingcode = mfilename('fullpath');

new_stimdrivenrate = exp( pstar(2)*log(stimdrivenrate.test) + pstar(1) );
PS                 = linearfilters.PostSpike.Filter;
xvalperformance = subR_xvalperformance_GLM(new_stimdrivenrate, testspikes_logicalraster, t_bin,PS)

fittedGLM.xvalperformance  = xvalperformance; 
eval(sprintf('save %s/%s.mat fittedGLM',glm_cellinfo.d_save,glm_cellinfo.cell_savename));

error('need new plotting  20015-07-02')
end


function NL_xvalperformance     = subR_xvalperformance_GLM( stimdrivenrate, logicalspike, t_bin,PS)
% AKHEITMAN 2015-06-24  it works!
% New version includes PS filter

params.bindur     = t_bin;
params.bins       = length(stimdrivenrate);
params.trials     = size(logicalspike,1);
params.testdur_seconds = params.bindur * params.bins ;   

% Set log-conditional as stim driven only
lcif_teststim = log(stimdrivenrate);
lcif = repmat(lcif_teststim, params.trials,1);  

lcif_ps = fastconv(logicalspike , [0; PS]', size(logicalspike,1), size(logicalspike,2) );    
lcif = lcif + lcif_ps;

glm_ratepersec  = exp(lcif);
glm_rateperbin  = params.bindur * glm_ratepersec;

spikerate_bin    = size(find(logicalspike(:))) /  size(logicalspike(:));      
model_null0      = spikerate_bin * ones(1, params.bins);
model_null       = repmat(model_null0, params.trials, 1);
null_logprob     = sum(eval_rasterlogprob(logicalspike, model_null, 'binary', 'conditioned'));
[raster_logprob_bin] = eval_rasterlogprob( logicalspike, glm_rateperbin,  'binary', 'conditioned') ;
glm_logprob       = sum(raster_logprob_bin);
glm_bits          = glm_logprob - null_logprob;
glm_bits_perspike = glm_bits / (sum(model_null0));
glm_bits_perbin   = glm_bits / params.bins;
glm_bits_persecond   = glm_bits / params.testdur_seconds;

NL_xvalperformance.note = 'Scores include optimized Non-Linearity';
NL_xvalperformance.logprob_null_raw            = null_logprob;
NL_xvalperformance.logprob_glm_raw      =  glm_logprob;
NL_xvalperformance.logprob_glm_bpspike  =  glm_bits_perspike;
NL_xvalperformance.logprob_glm_bpsec    =  glm_bits_persecond;

lcif_const  = lcif(1,:);
logical_sim = zeros(params.trials, params.bins);
% PS Filter 
cif_psgain = exp(PS);
ps_bins     = length(cif_psgain);
for i_trial = 1 : size(logicalspike,1)
    cif0         = exp(lcif_const);         
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

% IF NO PS FILTER
%{
for i_trial = 1 : size(logicalspike,1)
    cif         = exp(lcif_const);         
    binary_simulation = zeros(1,params.bins);
    for i = 1 : params.bins;
        roll = rand(1);
        if roll >  exp(-params.bindur*cif(i));
            binary_simulation(i)= 1;
        end
    end
    logical_sim(i_trial,:) = binary_simulation ;
end
%}
NL_xvalperformance.rasters.note           = 'glmsim includes altered non-linearity';
NL_xvalperformance.rasters.recorded       = logicalspike;
NL_xvalperformance.rasters.glm_sim        = logical_sim;
NL_xvalperformance.rasters.bintime        = params.bindur;
end   