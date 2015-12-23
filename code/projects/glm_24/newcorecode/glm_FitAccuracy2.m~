function glm_FitAccuracy2( Basepars, fit_type, Error_Only, simspertrial, movieoffset_fromzeroone )
% Version 2 2014-01-15
% Get rid of directory designed to include saving etc.
% Maybe just make this entirely obsolete?

%%% AK Heitman  2013-01-08

%%%  Takes 1-3 minutes depending on spike count   Only 10 millisec res.
%%% 1) simtorastermetrics_streamlined.m
%%%       Run Simulations, load or run the Variability of Ret Rasters
%%% 2) analyze_simtorastermetrics_streamlined.m
%%%      Compute Normalized SPKD between Raster 



[StimulusPars DirPars datarun_slv datarun_mas] = Directories_Params_v19_split(Basepars.exp_nm, fit_type, Basepars.map_type)

GLMPars  = Basepars.GLMPars;
fit_type = Basepars.GLMPars.fit_type; Head_CTYPE = Basepars.celltype;
d_save   = Basepars.d_save;

%%%%%% Settting Up Directories %%%%%%%%
plotRaster_dir         = sprintf('%s/Rasters',d_save);
AccuracyMetrics_dir = sprintf('%s/AccuracyMetrics',d_save);  clear d_sfx d_save fit_blks fit_type Head_CTYPE
if (~isdir(plotRaster_dir)),mkdir(plotRaster_dir); end
if (~isdir(AccuracyMetrics_dir)),mkdir(AccuracyMetrics_dir); end
%cd(AccuracyMetrics_dir);  

%%
%%%%%%%%%%%%%% RUN SPKD METRICS AND ANALYZE THE METRICS%%%%

[SimandMetrics] = simtorastermetrics_streamlined(Basepars,StimulusPars,DirPars,datarun_all,plotRaster_dir, computespikemetrics, simspertrial) %%% still a big long program

if ~SimandMetrics.MetricPars.debug, save(sprintf('%s/%s_SimulationsandMetrics.mat',AccuracyMetrics_dir, Basepars.fn_save), 'SimandMetrics');end
if  SimandMetrics.MetricPars.debug, save(sprintf('%s/%s_debugSimulationsandMetrics.mat',AccuracyMetrics_dir, Basepars.fn_save), 'SimandMetrics');end
if computespikemetrics.BW || computespikemetrics.NSEM    
    [MetricAnal]      = analyze_simtorastermetrics_streamlined(Basepars,SimandMetrics,computespikemetrics); %answer is just a dummy output  
   
    %%%%%%%%%%%% SAVE IN THE STANDARD OUTPUT FOLDER, AS WELL AS FOLDER THAT%%%
    %%%%%%%%%%%% COMPARES DIFFERENT PREPS %%%%%%%%
	if ~SimandMetrics.MetricPars.debug,save(sprintf('%s/%d_MetricAnal_%s.mat',AccuracyMetrics_dir, Basepars.headid,Basepars.fn_save), 'MetricAnal');end
    if SimandMetrics.MetricPars.debug, save(sprintf('%s/%d_debugMetricAnal_%s.mat',AccuracyMetrics_dir, Basepars.headid,Basepars.fn_save), 'MetricAnal');end
    Basepars_brief = Basepars;
    Basepars_brief = rmfield(Basepars_brief, 'trackprog_dir');%Basepars_brief = rmfield(Basepars_brief, 'kx_STA');
      %  Basepars_brief = rmfield(Basepars_brief, 'crop_idx');
    Basepars_brief = rmfield(Basepars_brief, 'p_evolve');
    clear SolSpace_MetricAnal; 
    SolSpace_MetricAnal_thisparset.debug_mode    = Basepars.debug_mode;  
    SolSpace_MetricAnal_thisparset.Basepars      = Basepars_brief;      SolSpace_MetricAnal_thisparset.MetricPars    = SimandMetrics.MetricPars;
    SolSpace_MetricAnal_thisparset.BW_Test       = MetricAnal.BW;      
    if isfield( MetricAnal, 'NSEM'),  SolSpace_MetricAnal_thisparset.NSEM_Test     = MetricAnal.NSEM; end
    SolSpace_MetricAnal_thisparset.fit_type      = Basepars.fit_type;
    SolSpace_MetricAnal_thisparset.cid           = Basepars.headid;     SolSpace_MetricAnal_thisparset.k_filtermode  = Basepars.k_filtermode;
    SolSpace_MetricAnal_thisparset.Coupling      = Basepars.Coupling;   SolSpace_MetricAnal_thisparset.BiDirect_CP   = Basepars.BiDirect_CP;
    SolSpace_MetricAnal_thisparset.tolfun        = GLMPars.tolfun;      SolSpace_MetricAnal_thisparset.binning       = GLMPars.binning;
    SolSpace_MetricAnal_thisparset.SolutionSpace = Basepars.SolutionSpace;
    if strcmp(Basepars.celltype,'ON-Parasol'),Ctype = 'ONPar';      end
    if strcmp(Basepars.celltype, 'OFF-Parasol'),Ctype = 'OFFPar';   end

    filename = sprintf('%s_%d_Metrics_ErrorOnly', Ctype ,  Basepars.headid);
    if isfield(computespikemetrics, 'Error_Only') && ~computespikemetrics.Error_Only
        filename = sprintf('%s_%d_Metrics_withintraSimvar', Ctype ,  Basepars.headid);
    end

    if  exist(sprintf('%s/%s.mat',FinalSpace_AnalysisDir,filename))
            load(sprintf('%s/%s.mat',FinalSpace_AnalysisDir,filename));  
            index = max(size(SolSpace_MetricAnal)) + 1;
            SolSpace_MetricAnal{index}  = SolSpace_MetricAnal_thisparset;
            clear index
    end
    if ~exist(sprintf('%s/%s.mat',FinalSpace_AnalysisDir,filename))
            SolSpace_MetricAnal{1}  = SolSpace_MetricAnal_thisparset;        
    end
    save(sprintf('%s/%s.mat',FinalSpace_AnalysisDir,filename), 'SolSpace_MetricAnal');   clear SolSpace_MetricAnal; 
end






%% COMPUTE THE LOG PROB
%keyboard
%%%% LOG PROB %%%%%
display('Computing Raster LogProb wrt Fitted GLM')
p_opt = Basepars.p_opt;
[Raster_LogProb]     = GLMcompute_rasterLogProbwrapper_streamlined(Basepars, SimandMetrics,StimulusPars, DirPars, datarun_all, p_opt);
if ~SimandMetrics.MetricPars.debug,save(sprintf('%s/%d_Raster_LogProb_%s.mat',AccuracyMetrics_dir, Basepars.headid,Basepars.fn_save), 'Raster_LogProb');end
if SimandMetrics.MetricPars.debug, save(sprintf('%s/%d_debugRaster_LogProb_%s.mat',AccuracyMetrics_dir, Basepars.headid,Basepars.fn_save), 'Raster_LogProb');end


if ~exist('Basepars_brief'), Basepars_brief = Basepars; end


SolSpace_LogProbAnal_thisparset.debug_mode    = Basepars.debug_mode;  
SolSpace_LogProbAnal_thisparset.Basepars = Basepars_brief;
SolSpace_LogProbAnal_thisparset.BW_Test  = Raster_LogProb.BW;
if isfield( Raster_LogProb, 'NSEM'), SolSpace_LogProbAnal_thisparset.NSEM_Test= Raster_LogProb.NSEM; end
SolSpace_LogProbAnal_thisparset.fit_type      = Basepars.fit_type;
SolSpace_LogProbAnal_thisparset.SolutionSpace = Basepars.SolutionSpace;
SolSpace_LogProbAnal_thisparset.Basepars      = Basepars_brief;      SolSpace_LogProbAnal_thisparset.MetricPars    = SimandMetrics.MetricPars;
SolSpace_LogProbAnal_thisparset.cid           = Basepars.headid;     SolSpace_LogProbAnal_thisparset.k_filtermode  = Basepars.k_filtermode;
SolSpace_LogProbAnal_thisparset.Coupling      = Basepars.Coupling;   SolSpace_LogProbAnal_thisparset.BiDirect_CP   = Basepars.BiDirect_CP;
SolSpace_LogProbAnal_thisparset.tolfun        = GLMPars.tolfun;      SolSpace_LogProbAnal_thisparset.binning       = GLMPars.binning;
clear SolSpace_LogProbAnal; 
%keyboard
if strcmp(Basepars.celltype,'ON-Parasol'),Ctype = 'ONPar';      end
if strcmp(Basepars.celltype, 'OFF-Parasol'),Ctype = 'OFFPar';   end
logprob_filename = sprintf('%s_%d_LogProb', Ctype ,  Basepars.headid);

%if ~Basepars.debug_mode
%{
if  exist(sprintf('%s/%s.mat',FinalSpace_AnalysisDir,logprob_filename))
        load(sprintf('%s/%s.mat',FinalSpace_AnalysisDir,logprob_filename));  
        index = max(size(SolSpace_LogProbAnal)) + 1;
        SolSpace_LogProbAnal{index}  = SolSpace_LogProbAnal_thisparset;
        clear index
end
if ~exist(sprintf('%s/%s.mat',FinalSpace_AnalysisDir,logprob_filename))
        SolSpace_LogProbAnal{1}  = SolSpace_LogProbAnal_thisparset;        
end
save(sprintf('%s/%s.mat',FinalSpace_AnalysisDir,logprob_filename), 'SolSpace_LogProbAnal');   clear SolSpace_LogProbAnal;
   %} 
    
end