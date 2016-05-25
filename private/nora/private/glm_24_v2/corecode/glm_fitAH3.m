function [p_opt,Basepars,f_opt] = glm_fitAH3(Basepars,spikesconcat_Home,spikesconcat_Neighbors,gopts,novelROIstim_Concat)
%%% Currently Only does Metric calculations on BW
%%% Major reconstructing on 2012-10-05
%%% IN SOME SENSE NOT MODULAR.. LOTS GET SAVED INTO VARIOUS FILES AND

%%% Integrate GLMPars

% isolated to a single neurons, it's spikes , outputdir, movie over ROI
% fminunc options, and some base params

% called by glm_AH_5

% This is a rip-off version of chaitu's mutiple_neuron.m.  Here the code is
% simplified to fit parameters of a *single* neuron (i.e., no coupling term
% between neurons).
% edoi@salk.edu, 2012-01-16,2012-04-16.

%calls create_histbasis
%calls pull_spikes
%calls grad_basis
%calls train_glm2 (rk2)
%calls nonsep_lgrad(raw)   train_glm_nonsep('raw'(

% new work and variable renaming for now

% NOT AN ACTUAL FUNCTION .. JUST CLEARS SPACES



keyboard
%%% stim and train pars need to be set   %%%
Stimpars.dt = Basepars.tstim; % this is the time (in sec) each stimulus frame is presented for.
Trainpars.baseneuron_idx = 1; % index of the neuron(s) being trained simultaneously
Trainpars.dt = Stimpars.dt / Basepars.spikebins_perstimframe; % time resolution of spikes    
Trainpars.neighbors_id = Basepars.cp_Neighbors;

neighbors_id = Basepars.cp_Neighbors;
% TOTAL MICROBINS  AND OFFSET %
microBins = Basepars.spikebins_perstimframe*Basepars.maxt;  % this is in (1/50)th units per frame
duration = Stimpars.dt * Basepars.maxt;
if (~isfield(Basepars,'frame_offset'))
   stim_offset = 0;
else
   stim_offset = Stimpars.dt*Basepars.frame_offset;
end
display('%%% Loading Spikes into Logical on Microbin Scale %%%')
%%% PULL SPIKES   %%%%%%%

Trainpars.logicalspike_microbin_Home      = sparse(logical(false(microBins,1)));         %in some sense initialize b4 puling spikes
Trainpars.logicalspike_microbin_Neighbors = sparse(logical(false(microBins,length(neighbors_id)))); 
%%% sp_times  is not really used,  but does give time in seconds for the
%%% conatenated blocks of interest
%%% Trainpars.D  has the action, in the fine time scale logical 
[Trainpars.sptimes_Home     ,Trainpars.logicalspike_microbin_Home     ,Trainpars.negSpikes_Home     ] = pull_spikesAH(spikesconcat_Home     ,Basepars.headid         ,microBins,duration,Trainpars.dt,stim_offset);
if Basepars.Coupling && length(Basepars.cp_Neighbors) >0
    [Trainpars.sptimes_Neighbors,Trainpars.logicalspike_microbin_Neighbors,Trainpars.negSpikes_Neighbors] = pull_spikesAH(spikesconcat_Neighbors,neighbors_id,microBins,duration,Trainpars.dt,stim_offset);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the postspike basis gradients (convolve each spike train with the basis fns)
%fprintf('Computing the spike train convolutions...');
%fprintf('done.\n');



%% purely a convlution of spikes with ps fileter   un optimized
%%% 30 secs for 7 cells.. this prolly doens't need to be optimized  yet 
display('%%% Convolving spikes with PS and CP basis vectors %%%%%%');
microBins_offset = Basepars.frame_offset*Basepars.spikebins_perstimframe;
Trainpars.psbasisGrad = grad_basisAH([Trainpars.negSpikes_Home; Trainpars.logicalspike_microbin_Home],Basepars.ps_basis,microBins_offset);
Trainpars.cpbasisGrad = [];

%keyboard
%%%% shift the 
if Basepars.Coupling && Basepars.BiDirect_CP
    oldNeighspike = Trainpars.logicalspike_microbin_Neighbors;
    couplingoffset = Basepars.cp_timebins;
    newCPpart1    = oldNeighspike( (couplingoffset+1):end , : );
    newCPpart2    = false( couplingoffset , size(newCPpart1,2) );
    newCP         = [newCPpart1 ; newCPpart2 ];
    
    Trainpars.logicalspike_microbin_Neighbors_preshift = oldNeighspike ; 
    Trainpars.logicalspike_microbin_Neighbors = newCP;
   
end
    


if Basepars.Coupling && length(Basepars.cp_Neighbors >0)
    Trainpars.cpbasisGrad = grad_basisAH([Trainpars.negSpikes_Neighbors; Trainpars.logicalspike_microbin_Neighbors],Basepars.cp_basis,microBins_offset);% note: coupling is not examined here.
end

%Nice and concatedneate   stim,   taking account for any offset
frame_idx = Basepars.frame_offset+1:Basepars.frame_offset+Basepars.maxt;   %%% total frames.. again   
Stimpars.movie_ROI = novelROIstim_Concat(:,frame_idx); clear novelROIstim_Concat  %%% the movie over the ROI

%%%%%%%%%%%% FIXED STA  or FIXED PS FILTER   %%%%%%%%%%%%%%%%%
if strcmp(Basepars.k_filtermode, 'STA')
    Basepars.kx_STA          = sum(fastconv(Stimpars.movie_ROI,Basepars.STA,Basepars.k_spacepixels,Basepars.maxt,Basepars.padval),1)';
end
if Basepars.ps_FIX
    ps_weights               = Basepars.ps_Filter;
    PS                       = Basepars.ps_basis*ps_weights; % Basepars.Mhist x N    
    Basepars.ps_FIXconvolved = (ps_weights') *cell2mat(Trainpars.psbasisGrad);
end
    
    

%keyboard

%% initial guess %%%
p0 = Basepars.p0;
Basepars.crop_idx = (1:Basepars.k_spacepixels)';
timeoffset = Basepars.frame_offset*Stimpars.dt; % offset time in seconds
%fprintf('Using frames %d-%d spanning the time range %f-%f sec.\n',...
%   Basepars.frame_offset+1,Basepars.frame_offset+Basepars.maxt,timeoffset,timeoffset+duration);


display('%%% Starting fminunc %%%%');
  %-- Fitting
switch Basepars.k_filtermode
   case 'rk2'
     % [p_opt f_opt g_opt,H_opt,exitflag,output] = train_glm2_AH(p0,Basepars,Stimpars,Trainpars,gopts,1);   %%% all good till here
       pstar  = p0;
       [pstar fstar eflag output] = fminunc(@(p) ll_func2AH(p,Basepars,Stimpars,Trainpars),pstar,gopts);
       [f_opt g_opt H_opt lcifs] = ll_func2AH(pstar,Basepars,Stimpars,Trainpars);
       p_opt = pstar;
   case 'STA'
       pstar  = p0;
       [pstar fstar eflag output] = fminunc(@(p) ll_func2AH(p,Basepars,Stimpars,Trainpars,Basepars.kx_STA),pstar,gopts);
       [f_opt g_opt H_opt lcifs] = ll_func2AH(pstar,Basepars,Stimpars,Trainpars,Basepars.kx_STA);
       p_opt = pstar;
   case 'raw'  % raw = nonsep; pixel- and frame-based filter, i.e., no assumption on K.
      Trainpars.lgrad = nonsep_lgrad(Basepars,stimpars,Trainpars);
      [p_opt,f_opt,g_opt,H_opt,exitflag,output] = train_glm_nonsep(p0,Basepars,Stimpars,Trainpars,gopts);
end

%% save & print
opt_param.p  = p_opt;
opt_param.p0 = p0;
opt_param.f  = f_opt;
opt_param.g  = g_opt;
opt_param.H  = H_opt;
opt_param.exitflag = eflag;
opt_param.output   = output;
opt_param.lcifs = lcifs;
Basepars.p_opt  = p_opt;
%Stimpars = rmfield(Stimpars,'movie_ROI');
Trainpars = rmfield(Trainpars,'psbasisGrad');
if(isfield(Trainpars,'lgrad'))
   Trainpars = rmfield(Trainpars,'lgrad');
end
if(isfield(Trainpars,'cpbasisGrad'))
   Trainpars = rmfield(Trainpars,'cpbasisGrad');
end
cid = Basepars.headid;
%final outputs into the folder fit

load(sprintf('%s/Track_Progress.mat',Basepars.trackprog_dir));
cd(Basepars.d_save);
save(fn_save,'Basepars','opt_param','Trainpars','cid','eflag','output','Track_Progress');
%keyboard;
print_glm_fitAH(cid,opt_param,Basepars,Trainpars,fn_save,Track_Progress);
%raster_glmAH;



%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% METRIC SPKD DIRECTORY SETTING %%%%  MAYBE THIS SHOULD GET WRAPPED AND
%%% MADE INTO ITS OWN FILE .. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[StimulusPars DirPars]  = Directories_Params_func_all(Basepars.exp_nm, Basepars.GLMPars.debug); 
GLMPars = Basepars.GLMPars;
fit_type = Basepars.GLMPars.fit_type; Head_CTYPE = Basepars.celltype;

d_save = Basepars.d_save;
if Basepars.Coupling && Basepars.BiDirect_CP, d_save = sprintf('%s/%s_FIT/%s/%s_ps%d_%s/%s',DirPars.output_dir,fit_type,Head_CTYPE,Basepars.k_filtermode,Basepars.ps_filternumber,'cpBiD',d_sfx);end
plotRaster_dir         = sprintf('%s/Rasters',d_save);
MetricAnalysis_dir = sprintf('%s/Metric_Analysis',d_save);  clear d_sfx d_save fit_blks fit_type Head_CTYPE
if (~isdir(plotRaster_dir)),mkdir(plotRaster_dir); end
if (~isdir(MetricAnalysis_dir)),mkdir(MetricAnalysis_dir); end
cd(MetricAnalysis_dir);  
SolutionSpace_AnalysisDir =sprintf('%s/SolutionSpace_Analysis',DirPars.output_dir);
if (~isdir(SolutionSpace_AnalysisDir)),mkdir(SolutionSpace_AnalysisDir); end    


%%%%%%%%%%%%%% RUN SPKD METRICS AND LOGPROB COMPUTATION%%%%
%keyboard
SimandMetrics = compute_simtoraster_metrics2(Basepars,opt_param,plotRaster_dir);
if ~SimandMetrics.MetricPars.debug, save(sprintf('%s/%d_SimulationsandMetrics_%s.mat',MetricAnalysis_dir, Basepars.headid,Basepars.fn_save), 'SimandMetrics');end
if SimandMetrics.MetricPars.debug,save(sprintf('%s/%d_debugSimulationsandMetrics_%s.mat',MetricAnalysis_dir, Basepars.headid,Basepars.fn_save), 'SimandMetrics');end

[MetricAnal]      = analyze_simtoraster_metrics2(Basepars,SimandMetrics); %answer is just a dummy output
if ~SimandMetrics.MetricPars.debug,save(sprintf('%s/%d_MetricAnal_%s.mat',MetricAnalysis_dir, Basepars.headid,Basepars.fn_save), 'MetricAnal');end
if SimandMetrics.MetricPars.debug, save(sprintf('%s/%d_dbugMetricAnal_%s.mat',MetricAnalysis_dir, Basepars.headid,Basepars.fn_save), 'MetricAnal');end
Basepars_brief = Basepars;
Basepars_brief = rmfield(Basepars_brief, 'trackprog_dir');Basepars_brief = rmfield(Basepars_brief, 'kx_STA');
Basepars_brief = rmfield(Basepars_brief, 'crop_idx');Basepars_brief = rmfield(Basepars_brief, 'p_evolve');
clear SolSpace_MetricAnal; 
SolSpace_MetricAnal_thisparset.debug_mode    = Basepars.debug_mode;  
SolSpace_MetricAnal_thisparset.Basepars      = Basepars_brief;      SolSpace_MetricAnal_thisparset.MetricPars    = SimandMetrics.MetricPars;
SolSpace_MetricAnal_thisparset.BW_Test       = MetricAnal.BW;      
if isfield( MetricAnal, 'NSEM'),  SolSpace_MetricAnal_thisparset.NSEM_Test     = MetricAnal.NSEM; end
SolSpace_MetricAnal_thisparset.cid           = Basepars.headid;     SolSpace_MetricAnal_thisparset.k_filtermode  = Basepars.k_filtermode;
SolSpace_MetricAnal_thisparset.Coupling      = Basepars.Coupling;   SolSpace_MetricAnal_thisparset.BiDirect_CP   = Basepars.BiDirect_CP;
SolSpace_MetricAnal_thisparset.tolfun        = GLMPars.tolfun;      SolSpace_MetricAnal_thisparset.binning       = GLMPars.binning;
SolSpace_MetricAnal_thisparset.SolutionSpace = Basepars.SolutionSpace;

if strcmp(Basepars.celltype,'ON-Parasol'),Ctype = 'ONPar';      end
if strcmp(Basepars.celltype, 'OFF-Parasol'),Ctype = 'OFFPar';   end
filename = sprintf('%s_%d_Metrics', Ctype ,  Basepars.headid);
if  exist(sprintf('%s/%s.mat',SolutionSpace_AnalysisDir,filename))
    load(sprintf('%s/%s.mat',SolutionSpace_AnalysisDir,filename));  
    index = max(size(SolSpace_MetricAnal)) + 1;
    SolSpace_MetricAnal{index}  = SolSpace_MetricAnal_thisparset;
    clear index
end
if ~exist(sprintf('%s/%s.mat',SolutionSpace_AnalysisDir,filename))
    SolSpace_MetricAnal{1}  = SolSpace_MetricAnal_thisparset;        
end

%if ~Basepars.debug_mode
    save(sprintf('%s/%s.mat',SolutionSpace_AnalysisDir,filename), 'SolSpace_MetricAnal');   clear SolSpace_MetricAnal;
%end
%%%% LOG PROB %%%%%
display('Computing Raster LogProb wrt Fitted GLM')
[Raster_LogProb]     = GLMcompute_rasterLogProbwrapper(Basepars, SimandMetrics, opt_param);
SolSpace_LogProbAnal_thisparset.debug_mode    = Basepars.debug_mode;  
SolSpace_LogProbAnal_thisparset.Basepars = Basepars_brief;
SolSpace_LogProbAnal_thisparset.BW_Test  = Raster_LogProb.BW;
if isfield( Raster_LogProb, 'NSEM'), SolSpace_LogProbAnal_thisparset.NSEM_Test= Raster_LogProb.NSEM; end
SolSpace_LogProbAnal_thisparset.SolutionSpace = Basepars.SolutionSpace;
SolSpace_LogProbAnal_thisparset.Basepars      = Basepars_brief;      SolSpace_LogProbAnal_thisparset.MetricPars    = SimandMetrics.MetricPars;
SolSpace_LogProbAnal_thisparset.cid           = Basepars.headid;     SolSpace_LogProbAnal_thisparset.k_filtermode  = Basepars.k_filtermode;
SolSpace_LogProbAnal_thisparset.Coupling      = Basepars.Coupling;   SolSpace_LogProbAnal_thisparset.BiDirect_CP   = Basepars.BiDirect_CP;
SolSpace_LogProbAnal_thisparset.tolfun        = GLMPars.tolfun;      SolSpace_LogProbAnal_thisparset.binning       = GLMPars.binning;
clear SolSpace_LogProbAnal; 
logprob_filename = sprintf('%s_%d_LogProb', Ctype ,  Basepars.headid);

%if ~Basepars.debug_mode
    if  exist(sprintf('%s/%s.mat',SolutionSpace_AnalysisDir,logprob_filename))
        load(sprintf('%s/%s.mat',SolutionSpace_AnalysisDir,logprob_filename));  
        index = max(size(SolSpace_LogProbAnal)) + 1;
        SolSpace_LogProbAnal{index}  = SolSpace_LogProbAnal_thisparset;
        clear index
    end
    if ~exist(sprintf('%s/%s.mat',SolutionSpace_AnalysisDir,logprob_filename))
        SolSpace_LogProbAnal{1}  = SolSpace_LogProbAnal_thisparset;        
    end
    save(sprintf('%s/%s.mat',SolutionSpace_AnalysisDir,logprob_filename), 'SolSpace_LogProbAnal');   clear SolSpace_LogProbAnal;
%end

%%%this is important!
delete(sprintf('%s/Track_Progress.mat',Basepars.trackprog_dir));

end

