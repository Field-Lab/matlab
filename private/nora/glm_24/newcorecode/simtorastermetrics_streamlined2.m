

%%% COMPLETE AND DOCUMENTED ON 2012-12-07   AK Heitman
%%% This can be a lengthy program

% version 2 will try to go without having both stimuli
% version 2 perhaps will be called for directly from glm_AH_19;

% 1) Computes Variability of Retina Raster
% 2) Run Model Simulations
% 3) Print Rasters and Model Raster Output


% Streamlined version of compute_simtoraster_metrics4
% Streamlined version complete  2013-01-08
% Only look at CentiSec Resolution
% Only 1 Simulation per 
% Version 3: GLMPars.computespikemetrics.BW and NSEM   
%            control which metrics we compute
%            2012-12-17
% version 2: no saving..  more modular

% SimandMetrics should automatically get saved into the correct directory 

% Inputs: Basepars and opt_param   (which will evetually get engulfed by Baseparams)
% Optional: plotRaster_dir :   program will save images of Rasters to
% specified output directory
% Outputs: SimandMetrics  large structure
%  .BW{i_trail,i_sim}  (also a corresponding NSEM)
%    .binnedlogical
%    .finalcif          conditinoal intensity func, including ps filters 
%    .blocknum          block associated with the trial
%    .simnum            simnum for the given trial
%    .sptimes_secs      model sim timing of spikes
%    .metric_spikes     model simulation spikes truncated for the metric
%    .ret_allspikes     the timing of RGC raster spikes on the given trial     
%    .ret_metricspikes  RGC spikes truncated that are used for metric 
%    .metricdist        distance to the Retina spike train of the sim
%           .millisec   param value for moving spikes in spkd
%           .halfcentiesec  .centisec .halfdecisec .decisec 5 param values
%    .normed_metricdist normalized by retina variability .. 
%                       * the more relevant metric *
% .BW_RasterVariation  (also exists a corresponding NSEM version)
%    .avgmillisec       Avg distance between trials of the retina raster
%    .stdmillisec       STD of distance between trial pairs of raster
% .MetricPars  
%    .debug             Cut the raster blocks and sims per block if TRUE
%    .startlate         How many seconds to ignore for metric computation
%    .endearly          Seconds on back end to ignore for metric comp
%    .sims_pertrial     how many simulations were run for each simulation
%                       set to 20

% This version is extremely modular... relatively well wrtten

% NORMED = DIVIDED BY RETINA RASTER VARIATION (AVG. TRIAL PAIR SPKD)

% SPKD is the key program.. JD Victor   

%CALLS:
% Directories_Params_func_all
% load_rasters   modular computation .. I believe
% raster_varation_metrics   modular computation .. I beleive
% logintensityfunction_lineafilterandbias
% logintensityfunction_coupling

function [SimandMetrics] = simtorastermetrics_streamlined2(Basepars,StimulusPars,DirPars,datarun_all,plotRaster_dir, simspertrial, stim_vrf)

GLMPars            = Basepars.GLMPars;
datarun_BW{1}      = datarun_all.master;  datarun_BW{2}   = datarun_all.BW_slv;
datarun_NSEM{1}    = datarun_all.master; datarun_NSEM{2}  = datarun_all.NSEM_slv; 

MetricPars.debug     = Basepars.debug_mode;
MetricPars.startlate = 1;  %% seconds of data we throw away in front
MetricPars.endearly  = 1;  %% seconds of data thrown away at the back
MetricPars.rasterplot= false;
MetricPars.sims_pertrial = simspertrial;

SimandMetrics.MetricPars = MetricPars;
%%%%%%%%%% Directory Fixing %%%%%%%%%%%
%%%%%%%%%%%%%%%%
clear Simulations; 

%%
for count = 1:2  
    %% Choose Verification Dateset and load Rasters
    tic
    clear datarun Slv_StimPars stim_vrf
    
    if count == 1
        stim_vrf = 'BW';   Slv_StimPars = StimulusPars.BW;    datarun = datarun_BW;
    elseif count == 2
        stim_vrf = 'NSEM'; Slv_StimPars = StimulusPars.NSEM;  datarun = datarun_NSEM;
    end
    Cell_ID = Basepars.headid; expnm_string = Basepars.exp_nm; logical_debug = GLMPars.debug; rastertype_string = stim_vrf;
    bins_perstimframe =  Basepars.spikebins_perstimframe;  % large enough so atmost  one spike per bin
    plot_logical = MetricPars.rasterplot;
    [~, Raster_sptimeSecs ] = load_rasters(Cell_ID, expnm_string, rastertype_string, bins_perstimframe, logical_debug, plot_logical);
    clear bins_perstimframe debug rastertype_string expnm_string Cell_ID  logical_debug plot_logical
        
    %% Load or Compute the Variability of the Retina Raster itself
    
    
        



    %%  LOAD THE MOVIES AND Sikes
    %%%% OUTPUT OF THIS SECTION TESTMOVIE_ROI  %%%%%%%
    %%%%  LOADING TESTRUN %%%%%%%%%
    load(sprintf('%s/%s/%s/prep_STAandROI/testrun_id%d_%dsq',DirPars.output_dir,'BW',Basepars.celltype,Basepars.headid,GLMPars.K_slen));  
    %%%%  LOADING MOVIES  %%%%%%%%%
    if count == 1
        testmovie0   = testrun.ROI.mov{1};
        frames_testmovie = size(testmovie0, 3);
        testmovie_ROI    = reshape(testmovie0, Basepars.k_spacepixels , frames_testmovie);
        clear frames_testmovie testmovie0
        
    elseif count == 2
        movies_dir = '/snle/home/snl-e/glm/output_AH';
        K_slen = GLMPars.K_slen;
        eval(sprintf('load %s/NSEM_Movie/eyemov1.mat',movies_dir));
        testmovie0 = X_rescaled;
        clear X; clear Xr; clear X_rescaled;
        movie_width  = Slv_StimPars.Width;
        movie_height = Slv_StimPars.Height;
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
    clear testrun;
    


   

    %% MODEL INPUT BUILDING   turn this also into a single function
    
    %try to keep time is y-axis convention
    % Trial independent portion
    p_opt = Basepars.p_opt;
    spikebins_perframe = Basepars.spikebins_perstimframe;  Paramind = Basepars.paramind;
    if strcmp(Basepars.k_filtermode, 'STA')
        LinearFilter        = p_opt(Basepars.paramind.L)*Basepars.STA;
        [lcif_mu, lcif_kx]  = logintensityfunction_linearfilterandbias (p_opt , Paramind , testmovie_ROI, spikebins_perframe, LinearFilter); 
    end
    if strcmp(Basepars.k_filtermode, 'rk2')
        [lcif_mu, lcif_kx]  = logintensityfunction_linearfilterandbias (p_opt , Paramind , testmovie_ROI, spikebins_perframe);
    end
    
    clear movie; clear spikebins_perframe; clear Paramind; clear LinearFilter
    
    p_ps     = p_opt(Basepars.paramind.PS);
    psfilter = Basepars.ps_basis * p_ps; 
    ps_cifgain  = exp(psfilter);

    if Basepars.Coupling && ~isempty(Basepars.cp_Neighbors) 
        [lcif_cp, ~] = logintensityfunction_coupling_streamlined (p_opt, Basepars, Slv_StimPars, datarun);
    end
    
    CIFcomponents.lcif_mu = lcif_mu;
    CIFcomponents.lcif_kx = lcif_kx;
    if  Basepars.Coupling  && ~isempty(Basepars.cp_Neighbors) 
        CIFcomponents.lcif_cp = lcif_cp;
    end
    
    if isfield(Basepars.paramind, 'SR')
        mov       = testmovie_ROI;
        onefilter = ones(Basepars.k_spacepixels , 10); 
        m_coit =  fastconv(mov,onefilter,Basepars.k_spacepixels,size(mov,2),Basepars.padval);
        pos_parts = find(m_coit >= 0) ;
        neg_parts = find(m_coit < 0 ) ;
        if strcmp(Basepars.celltype, 'ON-Parasol'), rec_m_coit = m_coit; rec_m_coit(neg_parts) = 0; end     
        if strcmp(Basepars.celltype, 'OFF-Parasol'), rec_m_coit = m_coit; rec_m_coit(pos_parts) = 0; end
        
        K = Basepars.STA;
        [~,mx_idx] = max(abs(K(:)));
        [~,~,mx_fr] = ind2sub([15,15,30],mx_idx);
        Spatial_Filter = Basepars.STA(:,mx_fr);
        
        Basepars.rast_rect = rec_m_coit;
        Basepars.rast_rect_conv_spSTA=  sum(fastconv(rec_m_coit, Spatial_Filter, Basepars.k_spacepixels,size(mov,2),Basepars.padval),1)';
        
        lcif_rec0 = p_opt(Basepars.paramind.SR) * Basepars.rast_rect_conv_spSTA;
        lcif_rec  = reprows(lcif_rec0, Basepars.spikebins_perstimframe);
        
        CIFcomponents.lcif_rec = lcif_rec;
    end
    CIFcomponents.cif_psgain = ps_cifgain;   
    
    
    
    clear lcif_cp lcif_mu lcif_kx p_ps psfilter ps_cifgain logical_BiDirect lcif_perneighbor 
    display('ModelLoaded into CIFcomponents structure')

    %% Simulations
    
    if  Basepars.Coupling  && ~isempty(Basepars.cp_Neighbors) 
        SimPars.Coupling                 = true;
        SimPars.trials                   = length(Slv_StimPars.RasterBlocks);
        SimPars.simulations_pertrial     = MetricPars.sims_pertrial;
    else
        SimPars.Coupling = false;
        SimPars.trials   = 1;
        SimPars.simulations_pertrial = length(RasterBlocks);
    end
    
    SimPars.seconds   = Slv_StimPars.nsec_o;
    SimPars.bin_dt    = Basepars.tstim / Basepars.spikebins_perstimframe;
    SimPars.PostSpike = true;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%  RUN SIMULATIONS %%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if count == 1, display('%%%%%%%Computing BW Simulations%%%%%%%');   end
    if count == 2, display('%%%%%%%%Computing NSEM Simulations%%%%%%%'); end
    Simulations = simulate_model(SimPars , CIFcomponents);
    if count == 1,   BW_Simulations = Simulations; end
    if count == 2, NSEM_Simulations = Simulations; end
    clear SimPars CIFcomponents         
   
    
    %%  Evaluating Metric Distance   Core of the Program!!!
    %   Calls spkd  which is located .. somewhere on the server 
    %   spkd is code of Reich and Victor
    %   last portion of the loop.. no calls other than spkd
    %   Modifies Simulations
    rastvar = Raster_Variability.summary_stats;
    if computespmetric 
        pctdonepre = 0;
        for i_trial = 1 : length(Slv_StimPars.RasterBlocks)
            retspikes        = Raster_sptimeSecs{i_trial}.spikes;
            ret_ind1         = find(retspikes > MetricPars.startlate ) ;
            ret_ind2         = find(retspikes < (Slv_StimPars.nsec_o-MetricPars.endearly) );
            ret_index        = intersect(ret_ind1 , ret_ind2);
            ret_metricspikes = retspikes(ret_index);
            for i_Sim = 1 : MetricPars.sims_pertrial
                simspikes         = Simulations{i_trial, i_Sim}.sptimes_secs;
                sim_ind1          = find(simspikes > MetricPars.startlate  ) ;
                sim_ind2          = find(simspikes < (Slv_StimPars.nsec_o-MetricPars.endearly));
                sim_index         = intersect(sim_ind1 , sim_ind2);
                sim_metricspikes  = simspikes(sim_index );   %%% CHnaged to real time spike train  starting at front off rather than 0 !!  11-28
                Simulations{i_trial,  i_Sim}.ret_allspikes            = retspikes;
                Simulations{i_trial,  i_Sim}.ret_metricspikes         = ret_metricspikes;
                Simulations{i_trial , i_Sim}.metric_spikes            = sim_metricspikes; 
                Simulations{i_trial , i_Sim}.metricdist.centisec      = spkd ( ret_metricspikes , sim_metricspikes , 100 );
                Simulations{i_trial , i_Sim}.normed_metricdist.centisec      = Simulations{i_trial , i_Sim}.metricdist.centisec / rastvar.avgcentisec;
            end
            pctdone =100* i_trial / length(Slv_StimPars.RasterBlocks);
            if round(pctdone/10) > round(pctdonepre/10);
                if count ==1, display(sprintf('Percent BW SimtoRasterMetric_%d',round(pctdone))); end
                if count ==2, display(sprintf('Percent NSEM SimtoRasterMetric_%d',round(pctdone))); end
            end
            pctdonepre = pctdone;
        end
        clear pctdone;   clear pctdonepre;  
        clear i_Sim i_trial
        clear  simspikes sim_ind1 sim_ind2 ret_ind1 ret_ind2 ret_index
    end
    
    if count == 1, SimandMetrics.BW   = Simulations; SimandMetrics.BW_RasterVariation   = rastvar; end
    if count == 2, SimandMetrics.NSEM = Simulations; SimandMetrics.NSEM_RasterVariation = rastvar; end
    clear Simulations rastvar    
        
    %% Diagnostic Plot
    % needs to be improved!   but good enough
 %   keyboard
    figure;
    if exist(plotRaster_dir)
        testseconds     = Slv_StimPars.nsec_o;
        simbins         = floor(Basepars.spikebins_perstimframe *testseconds /Basepars.tstim );
        plotseconds     = 2.5;
        figs            = floor (testseconds / plotseconds) ;
        plotbins        = floor(Basepars.spikebins_perstimframe *plotseconds /Basepars.tstim );

        if count == 1,  Sim_Raster = BW_Simulations; end
        if count == 2,  Sim_Raster = NSEM_Simulations; end

        for i_fig  = 1:figs-1 
            clf; 
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
            for k = 1:raster_trials
                x = x0 ( find(Raster_sptimeSecs{k}.binned_logical(bin_start:bin_end)) );
                y = (k/raster_trials) * ones(1, length(x) );
                plot(x,y,'r.','markersize',MS);

               % y = (k/raster_trials)*homespikes(bin_start:bin_end,k);
               % plot(x0,y,'r.','markersize',MS);
                if k == 1, hold on, ylabel('RGC'); end
                if k == raster_trials, hold off;    end
            end
            subplot(212)
            for k = 1:raster_trials
                x = x0 ( find(Sim_Raster{k,1}.binnedlogical(bin_start:bin_end)) );
                y = (k/raster_trials) * ones(1, length(x) );
                %y = (k/raster_trials)*Simulations{k}.spikes(bin_start:bin_end);
                plot(x,y,'k.','markersize',MS);
                if k == 1, hold on, ylabel('GLM'); xlabel(stim_vrf);end
                if k == raster_trials, hold off;    end

            end
            orient tall
            eval(sprintf('print -dpdf %s/%d_%s_part%d_%s',plotRaster_dir,Basepars.headid,stim_vrf,i_fig,Basepars.fn_save));
            hold off
        end
        clear x y k testseconds simbins plotseconds figs plotbins Sim_Raster x0 bin_end bin_start i_fig sec_start sec_end
    end
 
end
end





    