%%% Version 2  .. more Modular !!! !!!
%%% Version 1 This works but is an organisational mess as of 11-28

%cd(Basepars.d_save);  % take this off if we don't need it 
GLMPars                       = Basepars.GLMPars;
[StimulusPars DirPars ~]      = Directories_Params_func(Basepars.exp_nm, 'BW', GLMPars.debug);


%%%%%%%%%% Directory Fixing %%%%%%%%%%%
fit_type = GLMPars.fit_type; Head_CTYPE = Basepars.celltype;
switch (fit_type )
   case 'BW'
      fit_blks =   StimulusPars.BW_FitBlocks;
   case 'NSEM'
      fit_blks = StimulusPars.NSEM_FitBlocks;
end
if GLMPars.debug
    d_sfx = sprintf( 'debug_bin%d_blk%d',GLMPars.binning,length(fit_blks) );
else
    d_sfx = sprintf( 'bin%d_blk%d',GLMPars.binning,length(fit_blks) );
end
if Basepars.Coupling
	d_save = sprintf('%s/%s_FIT/%s/%s_ps%d_%s/%s',DirPars.output_dir,fit_type,Head_CTYPE,Basepars.k_filtermode,Basepars.ps_filternumber,'cpON',d_sfx);
end
if ~Basepars.Coupling
	d_save = sprintf('%s/%s_FIT/%s/%s_ps%d_%s/%s',DirPars.output_dir,fit_type,Head_CTYPE,Basepars.k_filtermode,Basepars.ps_filternumber,'cpOFF',d_sfx);
end
if Basepars.Coupling && Basepars.BiDirect_CP
	d_save = sprintf('%s/%s_FIT/%s/%s_ps%d_%s/%s',DirPars.output_dir,fit_type,Head_CTYPE,Basepars.k_filtermode,Basepars.ps_filternumber,'cpBiD',d_sfx);
end
Raster_dir         = sprintf('%s/Rasters',d_save);
MetricAnalysis_dir = sprintf('%s/Metric_Analysis',d_save); 
clear d_sfx d_save fit_blks fit_type Head_CTYPE
if (~isdir(Raster_dir))
    mkdir(Raster_dir); 
end
if (~isdir(MetricAnalysis_dir))
    mkdir(MetricAnalysis_dir); 
end
cd(MetricAnalysis_dir);  

%%%%%%%% Make Second Dir for Solution Space Analysis %%%%%%%%%
SolutionSpace_AnalysisDir =sprintf('%s/SolutionSpace_Analysis',DirPars.output_dir);
if (~isdir(SolutionSpace_AnalysisDir))
    mkdir(SolutionSpace_AnalysisDir); 
end



%%%%%%  Parameters concerning the Metrics used %%%%%%%%
MetricPars.startlate = 1;
MetricPars.endearly  = 1;
if GLMPars.debug
    MetricPars.sims_pertrial = 3; 
end
if ~GLMPars.debug
    MetricPars.sims_pertrial = 3;
end
clear Simulations; clear datarun; clear StimulusPars;
%%
for count = 1:2  
    %% Choose Verification Dateset and load Rasters
    if count == 1
        stim_vrf = 'BW';
    elseif count == 2
        stim_vrf = 'NSEM';
    end
    clear datarun
    [StimulusPars DirPars datarun]      = Directories_Params_func(Basepars.exp_nm, stim_vrf, GLMPars.debug);
       
    Cell_ID = Basepars.headid; expnm_string = Basepars.exp_nm; logical_debug = GLMPars.debug; rastertype_string = stim_vrf;
    bins_perstimframe =  Basepars.spikebins_perstimframe;  % large enough so atmost  one spike per bin
    [~, Raster_sptimeSecs ] = load_rasters(Cell_ID, expnm_string, rastertype_string, bins_perstimframe, logical_debug);
    clear bins_perstimframe debug rastertype_string expnm_string Cell_ID  logical_debug
    

    
    
    
    %% homespikes evaluation of avg distance etc.
    output_dir        = DirPars.output_dir; 
    CTYPE             = Basepars.celltype; 
    metvar_dirname    = sprintf('%s/%s/%s/metric_variability' ,output_dir , stim_vrf, CTYPE); 
    if (~isdir(metvar_dirname))
        mkdir(metvar_dirname); 
    end
    clear CTYPE; clear outputdir;
    trials = max ( size ( Raster_sptimeSecs) );
    % EVALUATE EACH PAIR METRIC DIST
    if ~exist(sprintf('%s/%d_RetRast_PairDist.mat',metvar_dirname,Basepars.headid))        
        RasterCell = Raster_sptimeSecs; raster_time = StimulusPars.nsec_o; start_delay = frontoff; end_early = backoff;
        maxpairs = 4 * trials;
        [Raster_Variability] = raster_variation_metrics( RasterCell, raster_time, start_delay, end_early, maxpairs);
        Raster_Variability.cid       = Basepars.headid;
        Raster_Variability.exp_nm    = Basepars.exp_nm;
        Raster_Variability.celltype  = Basepars.celltype
        clear RasterCell; clear raster_time; clear start_delay; clear end_early; clear maxpairs;
        display('%%% Computing and Saving Raster Variability Stats %%%')
        save(sprintf('%s/%d_RetRast_PairDist.mat',metvar_dirname,Basepars.headid),'Raster_Variability');
    else
        load(sprintf('%s/%d_RetRast_PairDist.mat',metvar_dirname,Basepars.headid))
        display('%%% Loaded Previously Computed Raster Variability Stats %%%')
    end
    clear metvar_dirname
    
        



    %%  LOAD THE MOVIES AND TESTRUN 
    %%%%  LOADING TESTRUN %%%%%%%%%
    load(sprintf('%s/%s/%s/prep_STAandROI/testrun_id%d_%dsq',DirPars.output_dir,GLMPars.fit_type,Basepars.celltype,Basepars.headid,GLMPars.K_slen));  
    %%%%  LOADING MOVIES  %%%%%%%%%
    if count == 1
        testmovie_ROI    = testrun.ROI.mov{1};
        frames_testmovie = size(testmovie_ROI, 3);
        testmovie_ROI    = reshape(testmovie_ROI, Basepars.k_spacepixels , frames_testmovie);
        clear frames_testmovie
        
    elseif count == 2
        movies_dir = '/snle/home/snl-e/glm/output_AH';
        K_slen = GLMPars.K_slen;
        eval(sprintf('load %s/NSEM_Movie/eyemov1.mat',movies_dir));
        testmovie0 = X_rescaled;
        clear X; clear Xr; clear X_rescaled;
        movie_width  = StimulusPars.NSEM_Width;
        movie_height = StimulusPars.NSEM_Height;
        testmovie1 = reshape(testmovie0, size(testmovie0,1), movie_width * movie_height);
        testmovie2 = transpose(testmovie1);
        testmovie  = reshape(testmovie2, movie_width, movie_height, size(testmovie0,1));        
        %%%%%%%%%%%%%%
        %%% Do something again with the datarun  Make official date to date
        %%% converter !!! %%%%%%%  look in prep folder for formula
        ROI_Xind = testrun.ROI.ROI_x;   %%% either way need the testrun with defined ROI ..
        ROI_Yind = testrun.ROI.ROI_y;
        testmovie_ROI = testmovie( ROI_Xind , ROI_Yind , : );
        frames_testmovie = size(testmovie_ROI, 3);
        testmovie_ROI = reshape(testmovie_ROI, Basepars.k_spacepixels , frames_testmovie);
        clear testmovie0; clear testmovie1; clear testmovie2;
        clear K_slen ROI_Xind ROI_Yind movie_height movie_width frames_testmovie 
    end
    % fast convolution needs to have this in Stimpars form
    testmovie    = testmovie_ROI;   clear testmovie_ROI; clear testrun;
    


   

    %% MODEL INPUT BUILDING   turn this also into a single function
    
    %try to keep time is y-axis convention
    % Trial independent portion
    p_opt = opt_param.p;
    spikebins_perframe = Basepars.spikebins_perstimframe;  Paramind = Basepars.paramind;
    if strcmp(Basepars.k_filtermode, 'STA')
        LinearFilter        = Basepars.STA;
        [lcif_mu, lcif_kx]  = logintensityfunction_linearfilterandbias (p_opt , Paramind , testmovie, spikebins_perframe, LinearFilter); 
    end
    if strcmp(Basepars.k_filtermode, 'rk2')
        [lcif_mu, lcif_kx] = logintensityfunction_linearfilterandbias (p_opt , Paramind , testmovie, spikebins_perframe);
    end
    clear movie; clear spikebins_perframe; clear Paramind; clear LinearFilter
 
    p_ps     = p_opt(Basepars.paramind.PS);
    psfilter = Basepars.ps_basis * p_ps; 
    ps_cifgain  = exp(psfilter);

    if Basepars.Coupling && ~isempty(Basepars.cp_Neighbors) 
        spikebins_perblock = Basepars.spikebins_perstimframe * (StimulusPars.ntb_o-1) * StimulusPars.frames_pertrigger;
        block_duration     = StimulusPars.nsec_o;
      %  spikebins      = frames_testmovie * Basepars.spikebins_perstimframe;
      %  binrate_persec = Basepars.tstim / Basepars.spikebins_perstimframe;
        Paramind  = Basepars.paramind; Neighbor_IDs = Basepars.cp_Neighbors  ; %% should be all there.. all misses should be eliminated
        Coupling_Basis = Basepars.cp_basis;
        if count == 1
            Raster_Blocks = StimulusPars.BW_RasterBlocks;
        elseif count ==2
            Raster_Blocks = StimulusPars.NSEM_RasterBlocks;
        end  
        logical_BiDirect=Basepars.BiDirect_CP;
        [lcif_cp, ~] = logintensityfunction_coupling (p_opt, datarun, Paramind ,Neighbor_IDs, Coupling_Basis, Raster_Blocks, spikebins_perblock,block_duration, logical_BiDirect); 
    end
    
    CIFcomponents.lcif_mu = lcif_mu;
    CIFcomponents.lcif_kx = lcif_kx;
    if  Basepars.Coupling  && ~isempty(Basepars.cp_Neighbors) 
        CIFcomponents.lcif_cp = lcif_cp;
    end
    CIFcomponents.cif_psgain = ps_cifgain;   
    clear lcif_cp lcif_mu lcif_kx p_ps psfilter ps_cifgain logical_BiDirect
    display('ModelLoaded into CIFcomponents structure')

    %% Simulations
    
    if  Basepars.Coupling  && ~isempty(Basepars.cp_Neighbors) 
        SimPars.Coupling                 = true;
        SimPars.trials                   = length(Raster_Blocks);
        SimPars.simulations_pertrial     = MetricPars.sims_pertrial;
    else
        SimPars.Coupling = false;
        SimPars.trials   = 1;
        SimPars.simulations_pertrial = length(Raster_Blocks);
    end
    SimPars.seconds   = StimulusPars.nsec_o;
    SimPars.bin_dt    = Basepars.tstim / Basepars.spikebins_perstimframe;
    SimPars.PostSpike = true;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%  RUN SIMULATIONS %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Simulations = simulate_model(SimPars , CIFcomponents);
    clear SimPars CIFcomponents            
   
    
    %%  Evaluating Metric Distance   Core of the Program!!!
    pctdonepre = 0;
    for i_trial = 1 : length(Raster_Blocks)
        retspikes        = Raster_sptimeSecs{i_trial}.spikes;
        ret_ind1         = find(retspikes > MetricPars.startlate ) ;
        ret_ind2         = find(retspikes < (StimulusPars.nsec_o-MetricPars.endearly) );
        ret_index        = intersect(ret_ind1 , ret_ind2);
        ret_metricspikes = retspikes(ret_index);
        
        for i_Sim = 1 : MetricPars.sims_pertrial
            simspikes         = Simulations{i_trial, i_Sim}.sptimes_secs;
            sim_ind1          = find(simspikes > MetricPars.startlate  ) ;
            sim_ind2          = find(simspikes < (StimulusPars.nsec_o-MetricPars.endearly));
            sim_index         = intersect(sim_ind1 , sim_ind2);
            sim_metricspikes  = simspikes(sim_index );   %%% CHnaged to real time spike train  starting at front off rather than 0 !!  11-28
            
            Simulations{i_trial,  i_Sim}.ret_metricspikes         = ret_metricspikes;
            Simulations{i_trial , i_Sim}.metric_spikes            = sim_metricspikes; 
            Simulations{i_trial , i_Sim}.metricdist_millisec      = spkd ( ret_metricspikes , sim_metricspikes , 1000);
            Simulations{i_trial , i_Sim}.metricdist_halfcentisec  = spkd ( ret_metricspikes , sim_metricspikes , 200 );
            Simulations{i_trial , i_Sim}.metricdist_centisec      = spkd ( ret_metricspikes , sim_metricspikes , 100 );
            Simulations{i_trial , i_Sim}.metricdist_halfdecise    = spkd ( ret_metricspikes , sim_metricspikes , 20  );    
            Simulations{i_trial , i_Sim}.metricdist_decisec       = spkd ( ret_metricspikes , sim_metricspikes , 10  );    
        end
        
        pctdone =100* i_trial / trials;
        if round(pctdone/10) > round(pctdonepre/10);
            display(sprintf('PercentSimtoRasterMetric_%d',round(pctdone)));
        end
        pctdonepre = pctdone;
    end
    clear pctdone;   clear pctdonepre;  
    
    if count == 1
        BW_Simulations    = Simulations;
    elseif count == 2
        NSEM_Simulations  = Simulations;
    end
    clear Simulations; clear datarun; clear StimulusPars;
end
save(sprintf('%s/%d_SimulationsandMetrics_%s.mat',MetricAnalysis_dir, Basepars.headid,Basepars.fn_save), 'BW_Simulations' , 'NSEM_Simulations');%
  % eval(sprintf(
%{     
%% Analyze Metrics   make some figures   
     
     % per trial distance 
     clear MetricAnalysis
     MetricAnalyis.pertrial    = cell (trials , 1);
     
     for i_trial = 1 : trials
        atemp = [Simulations{ i_trial , :}]; 
        MetricAnalysis.pertrial{i_trial}.avgdist_millisec       =  mean ([atemp.metricdist_millisec]); 
        MetricAnalysis.pertrial{i_trial}.stddist_millisec       =  std ([atemp.metricdist_millisec]); 
        MetricAnalysis.pertrial{i_trial}.norm_avgdist_millisec  =  mean ([atemp.metricdist_millisec]) / ret_metricvar.avgmillisec; 
        MetricAnalysis.pertrial{i_trial}.norm_stddist_millisec  =  std ([atemp.metricdist_millisec]) / (ret_metricvar.avgmillisec); 
        
        MetricAnalysis.pertrial{i_trial}.avgdist_mcsec          =  mean ([atemp.metricdist_mcsec]); 
        MetricAnalysis.pertrial{i_trial}.stddist_mcsec          =  std ([atemp.metricdist_mcsec]); 
        MetricAnalysis.pertrial{i_trial}.norm_avgdist_mcsec     =  mean ([atemp.metricdist_mcsec]) / ret_metricvar.avgmcsec; 
        MetricAnalysis.pertrial{i_trial}.norm_stddist_mcsec     =  std ([atemp.metricdist_mcsec]) / (ret_metricvar.avgmcsec); 
        
        MetricAnalysis.pertrial{i_trial}.avgdist_centisec       =  mean ([atemp.metricdist_centisec]); 
        MetricAnalysis.pertrial{i_trial}.stddist_centisec       =  std ([atemp.metricdist_centisec]);
        MetricAnalysis.pertrial{i_trial}.norm_avgdist_centisec  =  mean ([atemp.metricdist_centisec]) / ret_metricvar.avgcentisec; 
        MetricAnalysis.pertrial{i_trial}.norm_stddist_centisec  =  std ([atemp.metricdist_centisec]) / (ret_metricvar.avgcentisec); 
        
        MetricAnalysis.pertrial{i_trial}.avgdist_cdsec          =  mean ([atemp.metricdist_cdsec]); 
        MetricAnalysis.pertrial{i_trial}.stddist_cdsec          =  std ([atemp.metricdist_cdsec]); 
        MetricAnalysis.pertrial{i_trial}.norm_avgdist_cdsec     =  mean ([atemp.metricdist_cdsec]) / ret_metricvar.avgcdsec; 
        MetricAnalysis.pertrial{i_trial}.norm_stddist_cdsec     =  std ([atemp.metricdist_cdsec]) / (ret_metricvar.avgcdsec); 
        
        MetricAnalysis.pertrial{i_trial}.avgdist_decisec        =  mean ([atemp.metricdist_decisec]); 
        MetricAnalysis.pertrial{i_trial}.stddist_decisec        =  std ([atemp.metricdist_decisec]);
        MetricAnalysis.pertrial{i_trial}.norm_avgdist_decisec   =  mean ([atemp.metricdist_decisec]) / ret_metricvar.avgdecisec; 
        MetricAnalysis.pertrial{i_trial}.norm_stddist_decisec   =  std ([atemp.metricdist_decisec]) / (ret_metricvar.avgdecisec); 
     end
     MetricAnalysis.retinavariability = ret_metricvar;
     MetricAnalysis.cellid = cid;   
     
     atemp = [MetricAnalysis.pertrial{:}];
         
     
     MetricAnalysis.crosstrial.avgdistfromret_norm_milli = mean( [atemp.norm_avgdist_millisec] ); 
     MetricAnalysis.crosstrial.avgdistfromret_norm_mc = mean( [atemp.norm_avgdist_mcsec] ); 
     MetricAnalysis.crosstrial.avgdistfromret_norm_centi = mean( [atemp.norm_avgdist_centisec] ); 
     MetricAnalysis.crosstrial.avgdistfromret_norm_cd = mean( [atemp.norm_avgdist_cdsec] ); 
     MetricAnalysis.crosstrial.avgdistfromret_norm_deci = mean( [atemp.norm_avgdist_decisec] ); 
  
     MetricAnalysis.crosstrial.avgstddistfromret_norm_milli = mean( [atemp.norm_stddist_millisec] ); 
     MetricAnalysis.crosstrial.avgstddistfromret_norm_mc = mean( [atemp.norm_stddist_mcsec] ); 
     MetricAnalysis.crosstrial.avgstddistfromret_norm_centi = mean( [atemp.norm_stddist_centisec] ); 
     MetricAnalysis.crosstrial.avgstddistfromret_norm_cd = mean( [atemp.norm_stddist_cdsec] ); 
     MetricAnalysis.crosstrial.avgstddistfromret_norm_deci = mean( [atemp.norm_stddist_decisec] ); 
     
     
     
     clear choseones chosenpairs
     
     simnum         = size(Simulations,2);
     simpairs       = nchoosek(1:simnum, 2);
     totalsimpairs  = min (20, size(simpairs,1));
     if GLMPars.debug
         chosenones = 1;
     end
     if ~GLMPars.debug
         chosenones     = ceil ( size(simpairs,1)* rand (totalsimpairs , 1) ) ;
     end
     chosenpairs    = simpairs(chosenones,:);
     
     clear chosenones
     MetricAnalysis.intratrial.simpairdist = cell(trials, 1);
     clear i_pair; clear i_trial
     pctdonepre =0;
     for i_trial = 1 : trials
         MetricAnalysis.intratrial.simpairdist{i_trial}.allsimpairs_norm_milli = zeros(totalsimpairs, 1);
         MetricAnalysis.intratrial.simpairdist{i_trial}.allsimpairs_norm_mc    = zeros(totalsimpairs, 1);
         MetricAnalysis.intratrial.simpairdist{i_trial}.allsimpairs_norm_centi = zeros(totalsimpairs, 1);
         MetricAnalysis.intratrial.simpairdist{i_trial}.allsimpairs_norm_cd    = zeros(totalsimpairs, 1);
         MetricAnalysis.intratrial.simpairdist{i_trial}.allsimpairs_norm_deci  = zeros(totalsimpairs, 1);
         for i_pair = 1: totalsimpairs
            train1 = Simulations{i_trial , chosenpairs(i_pair,1)}.spikesecs_frontbackoff;
            train2 = Simulations{i_trial , chosenpairs(i_pair,2)}.spikesecs_frontbackoff;           
            MetricAnalysis.intratrial.simpairdist{i_trial}.allsimpairs_norm_milli(i_pair)     = spkd ( train1 , train2 , 1000) / ret_metricvar.avgmillisec;
            MetricAnalysis.intratrial.simpairdist{i_trial}.allsimpairs_norm_mc(i_pair)        = spkd ( train1 , train2 , 500 ) / ret_metricvar.avgmcsec;
            MetricAnalysis.intratrial.simpairdist{i_trial}.allsimpairs_norm_centi(i_pair)     = spkd ( train1 , train2 , 100 ) / ret_metricvar.avgcentisec; 
            MetricAnalysis.intratrial.simpairdist{i_trial}.allsimpairs_norm_cd(i_pair)        = spkd ( train1 , train2 , 50  ) / ret_metricvar.avgcdsec;
            MetricAnalysis.intratrial.simpairdist{i_trial}.allsimpairs_norm_deci(i_pair)      = spkd ( train1 , train2 , 10  ) / ret_metricvar.avgdecisec;
         end
        pctdone =100* i_trial / trials;
        if round(pctdone/10) > round(pctdonepre/10);
            display(sprintf('PercentSimulationVaraibilityDone_%d',round(pctdone)));
        end
        pctdonepre = pctdone;
         
     end
     clear pctdonepre pctdone
     clear atemp
     atemp = [MetricAnalysis.intratrial.simpairdist{:}];
     
     MetricAnalysis.intratrial.avgsimpairdist_norm_milli = mean( [atemp.allsimpairs_norm_milli] );
     MetricAnalysis.intratrial.avgsimpairdist_norm_mc    = mean( [atemp.allsimpairs_norm_mc   ] ); 
     MetricAnalysis.intratrial.avgsimpairdist_norm_centi = mean( [atemp.allsimpairs_norm_centi] ); 
     MetricAnalysis.intratrial.avgsimpairdist_norm_cd    = mean( [atemp.allsimpairs_norm_cd   ] ); 
     MetricAnalysis.intratrial.avgsimpairdist_norm_deci  = mean( [atemp.allsimpairs_norm_deci ] ); 

     %%%%%%%%%%%%%%%%%
     %%%%%%%%%%%%%%%%%
     %MetricAnalysis.nsec = 
     
     display('madeit')
     
     

     
     
  
     
     if count ==1
         BW_MetricAnal  = MetricAnalysis;
     end
     if count ==2
         NSEM_MetricAnal = MetricAnalysis;
     end
    clear MetricAnalysis
    
    
    
    
    
    
    cd ../Rasters
    testseconds     = nsec_o;
    plotseconds     = 2.5;
    figs            = floor (testseconds / plotseconds) ;
    plotbins        = floor(Basepars.spikebins_perstimframe *plotseconds /Basepars.tstim );
    for i_fig  = 1 : figs
        figure, clf
        MS = 6;
        bin_start = 1 + (i_fig-1) * plotbins;
        bin_end   = i_fig * plotbins ;   
        sec_start = 0 + (i_fig-1) *plotseconds;
        sec_end   = i_fig *plotseconds;
        x0         = linspace( sec_start , sec_end , plotbins);

        if bin_end > simbins
            bin_end = simbins;
            x0 =  linspace( sec_start , sec_end ,bin_end - bin_start +1);
        end

        %%%%%%%%%%%%
        subplot(211);
        for k = 1:trials
            x = x0 ( find(homespikes(bin_start:bin_end,k)) );
            y = (k/trials) * ones(1, length(x) );
            plot(x,y,'r.','markersize',MS);

           % y = (k/trials)*homespikes(bin_start:bin_end,k);
           % plot(x0,y,'r.','markersize',MS);
            if k == 1, hold on, end
        end
        subplot(212)
        for k = 1:trials
            x = x0 ( find(Simulations{k,1}.spikes(bin_start:bin_end)) );
            y = (k/trials) * ones(1, length(x) );
            %y = (k/trials)*Simulations{k}.spikes(bin_start:bin_end);
            plot(x,y,'k.','markersize',MS);
            if k == 1, hold on, end
        end
        orient tall
        eval(sprintf('print -dpdf %d_%s_part%d_%s',Basepars.headid,stim_vrf,i_fig,Basepars.fn_save));
    end
    cd ../Metric_Analysis
    
    if count == 1
         BW_Simulations = Simulations;
    end
    if count == 2
        NSEM_Simulations = Simulations;
    end
    
    clear datarun; clear Simulations; clear Pairwise; clear ret_metricvar ; clear retsp_metricspikes; clear trials

end
save(sprintf('%d_MetricAnal_%s',cid,Basepars.fn_save), 'BW_MetricAnal' , 'NSEM_MetricAnal');% 'BW_Simulations','NSEM_Simulations');

if strcmp(Basepars.celltype,'ON-Parasol')
    Ctype = 'ONPar';
elseif strcmp(Basepars.celltype, 'OFF-Parasol')
    Ctype = 'OFFPar';
end
filename = sprintf('%s_%d_Metrics', Ctype ,  Basepars.headid)

if exist(sprintf('%s/%s.mat',MetricAnalysis_Dir,filename))
    load(sprintf('%s/%s.mat',MetricAnalysis_Dir,filename));
    
    index = max(size(MetricAnal)) + 1;
    
    MetricAnal{index}.cid           = Basepars.headid;
    MetricAnal{index}.BW_Test       = BW_MetricAnal;
    MetricAnal{index}.NSEM_Test     = NSEM_MetricAnal;
  %  MetricAnal{index}.BW_Sim        = BW_Simulations;
  %  MetricAnal{index}.NSEM_Sim      = NSEM_Simulations;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MetricAnal{index}.Coupling      = Basepars.Coupling;
    MetricAnal{index}.BiDirect_CP   = Basepars.BiDirect_CP;
    MetricAnal{index}.ps_FIX        = Basepars.ps_FIX;
    MetricAnal{index}.tolx          = GLMPars.tolx;
    MetricAnal{index}.tolfun        = GLMPars.tolfun;
    MetricAnal{index}.binning       = GLMPars.binning;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MetricAnal{index}.Basepars      = Basepars;

end

if ~exist(sprintf('%s/%s.mat',MetricAnalysis_Dir,filename))
    
    MetricAnal{1}.cid           = Basepars.headid;
    MetricAnal{1}.BW_Test       = BW_MetricAnal;
    MetricAnal{1}.NSEM_Test     = NSEM_MetricAnal;
  %  MetricAnal{1}.BW_Sim        = BW_Simulations;
  %  MetricAnal{1}.NSEM_Sim      = NSEM_Simulations;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MetricAnal{1}.Coupling      = Basepars.Coupling;
    MetricAnal{1}.BiDirect_CP   = Basepars.BiDirect_CP;
    MetricAnal{1}.ps_FIX        = Basepars.ps_FIX;
    MetricAnal{1}.tolx          = GLMPars.tolx;
    MetricAnal{1}.tolfun        = GLMPars.tolfun;
    MetricAnal{1}.binning       = GLMPars.binning;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MetricAnal{1}.Basepars      = Basepars
    

end


save(sprintf('%s/%s.mat',MetricAnalysis_Dir,filename), 'MetricAnal');
 



%}

%% Plots for Comparing BW and NSEM runs





    