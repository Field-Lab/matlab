clear; clear -global; close all
runtype = 'network'; %runtype = 'local';
 
% INITIALIZATION AND DIRECTORY IDENTIFICATION / IMPORTANT PARAMS
local_homedir ='/Users/akheitman/Matlab_code/glm_AH_current';
network_currentdir = '/netapp/snle/home/snl-e/matlab-standard/code/projects/call_glm_AH';
network_homedir ='/netapp/snle/home/snl-e/matlab-standard/code/projects/glm_AH_15b'; 
if strcmp(runtype,   'local'), homedir =   local_homedir; end
if strcmp(runtype, 'network'), homedir = network_homedir; end
path(homedir, path);
if strcmp(runtype,'network'), path(network_homedir,path); end

matfiledir = '/netapp/snle/home/snl-e/glm/output_AH/2012-08-09-3/BW_Fit/STA_ps36_cpBiD/bin20_blk52_tolfun5_cpbasisno8'
cids = [1471 1786 3676 5086 5161]; ctype ='OFFPar';
cells = length(cids);
filenames = cell(1,cells);
for i_cell = 1 : cells
    filenames{i_cell} = sprintf('%s/%s_%d.mat', matfiledir, ctype, cids(i_cell));
end

computespikemetrics.BW = true;
computespikemetrics.NSEM = true;
computespikemetrics.Error_Only = true;


for i_cell = 1: cells
    load( filenames{i_cell}, 'Basepars','opt_param');
    
    [StimulusPars DirPars]  = Directories_Params_func_all(Basepars.exp_nm, Basepars.GLMPars.debug); 
    GLMPars = Basepars.GLMPars;
    fit_type = Basepars.GLMPars.fit_type; Head_CTYPE = Basepars.celltype;

    d_save = Basepars.d_save;
    %if Basepars.Coupling && Basepars.BiDirect_CP, d_save = sprintf('%s/%s_FIT/%s/%s_ps%d_%s/%s',DirPars.output_dir,fit_type,Head_CTYPE,Basepars.k_filtermode,Basepars.ps_filternumber,'cpBiD',d_sfx);end
    plotRaster_dir         = sprintf('%s/Rasters',d_save);
    MetricAnalysis_dir = sprintf('%s/Metric_Analysis',d_save);  clear d_sfx d_save fit_blks fit_type Head_CTYPE
    if (~isdir(plotRaster_dir)),mkdir(plotRaster_dir); end
    if (~isdir(MetricAnalysis_dir)),mkdir(MetricAnalysis_dir); end
    cd(MetricAnalysis_dir);  
    SolutionSpace_AnalysisDir =sprintf('%s/SolutionSpace_Analysis',DirPars.output_dir);
    if Basepars.debug_mode, SolutionSpace_AnalysisDir =sprintf('%s/dbug_SolutionSpace_Analysis',DirPars.output_dir); end
    if (~isdir(SolutionSpace_AnalysisDir)),mkdir(SolutionSpace_AnalysisDir); end    


    %%%%%%%%%%%%%% RUN SPKD METRICS AND LOGPROB COMPUTATION%%%%
    [SimandMetrics] = compute_simtoraster_metrics3(Basepars,plotRaster_dir, computespikemetrics);
    if ~SimandMetrics.MetricPars.debug, save(sprintf('%s/%d_SimulationsandMetrics_%s.mat',MetricAnalysis_dir, Basepars.headid,Basepars.fn_save), 'SimandMetrics');end
    if  SimandMetrics.MetricPars.debug, save(sprintf('%s/%d_debugSimulationsandMetrics_%s.mat',MetricAnalysis_dir, Basepars.headid,Basepars.fn_save), 'SimandMetrics');end



    if computespikemetrics.BW || computespikemetrics.NSEM

        [MetricAnal]      = analyze_simtoraster_metrics3(Basepars,SimandMetrics,computespikemetrics); %answer is just a dummy output   
        if ~SimandMetrics.MetricPars.debug,save(sprintf('%s/%d_MetricAnal_%s.mat',MetricAnalysis_dir, Basepars.headid,Basepars.fn_save), 'MetricAnal');end
        if SimandMetrics.MetricPars.debug, save(sprintf('%s/%d_debugMetricAnal_%s.mat',MetricAnalysis_dir, Basepars.headid,Basepars.fn_save), 'MetricAnal');end
        Basepars_brief = Basepars;
        Basepars_brief = rmfield(Basepars_brief, 'trackprog_dir');%Basepars_brief = rmfield(Basepars_brief, 'kx_STA');
      %  Basepars_brief = rmfield(Basepars_brief, 'crop_idx');
        Basepars_brief = rmfield(Basepars_brief, 'p_evolve');
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

        filename = sprintf('%s_%d_Metrics_ErrorOnly', Ctype ,  Basepars.headid);
        if isfield(computespikemetrics, 'Error_Only') && ~computespikemetrics.Error_Only
            filename = sprintf('%s_%d_Metrics_withintraSimvar', Ctype ,  Basepars.headid);
        end

        if  exist(sprintf('%s/%s.mat',SolutionSpace_AnalysisDir,filename))
            load(sprintf('%s/%s.mat',SolutionSpace_AnalysisDir,filename));  
            index = max(size(SolSpace_MetricAnal)) + 1;
            SolSpace_MetricAnal{index}  = SolSpace_MetricAnal_thisparset;
            clear index
        end
        if ~exist(sprintf('%s/%s.mat',SolutionSpace_AnalysisDir,filename))
            SolSpace_MetricAnal{1}  = SolSpace_MetricAnal_thisparset;        
        end
        save(sprintf('%s/%s.mat',SolutionSpace_AnalysisDir,filename), 'SolSpace_MetricAnal');   clear SolSpace_MetricAnal;
    end


    %% COMPUTE THE LOG PROB

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
end