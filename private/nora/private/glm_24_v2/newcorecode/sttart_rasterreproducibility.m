exp_nm = '2013-08-09-3';
[StimulusPars DirPars doutput_dir        = DirPars.output_dir; 
    



CTYPE             = Basepars.celltype; 
    metvar_dirname    = sprintf('%s/%s/%s/metric_variability' ,output_dir , stim_vrf, CTYPE); 
    if (~isdir(metvar_dirname))
        mkdir(metvar_dirname); 
    end
    clear CTYPE; clear outputdir;
    raster_trials = max ( size ( Raster_sptimeSecs) );

    if ~exist(sprintf('%s/%d_RetRast_PairDist.mat',metvar_dirname,Basepars.headid))        
        RasterCell = Raster_sptimeSecs;  start_delay = frontoff; end_early = backoff;
        if strcmp(stim_vrf, 'BW')  , raster_time = Slv_StimPars.BW.nsec_o;   end
        if strcmp(stim_vrf, 'NSEM'), raster_time = Slv_StimPars.NSEM.nsec_o; end
        
        maxpairs = 4 * raster_trials;
        [Raster_Variability]         = raster_variation_metrics_streamlined( RasterCell, raster_time, start_delay, end_early, maxpairs);
        Raster_Variability.cid       = Basepars.headid;
        Raster_Variability.exp_nm    = Basepars.exp_nm;
        Raster_Variability.celltype  = Basepars.celltype;
        clear RasterCell; clear raster_time; clear start_delay; clear end_early; clear maxpairs;
        if count == 1, display('%%% Computing and Saving BW Raster Variability Stats %%%'); end
        if count == 2, display('%%% Computing and Saving NSEM Raster Variability Stats %%%'); end
        save(sprintf('%s/%d_RetRast_PairDist.mat',metvar_dirname,Basepars.headid),'Raster_Variability');
    else
        load(sprintf('%s/%d_RetRast_PairDist.mat',metvar_dirname,Basepars.headid))
        display('%%% Loaded Previously Computed Raster Variability Stats %%%')
    end
    clear metvar_dirname atarun_slv datarun_mas] = Directories_Params_v19_split(exp_nm, GLMPars.fit_type, map_type)