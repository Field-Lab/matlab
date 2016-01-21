% AK HEITMAN
% Reviewed 2012-12-13
%
% Wrapper that calls "compute_rasterLogProb"
% Specifically in a load output Basepars and p_opt  context
% Also calls Directoies_Params_func_all

%INPUT: Basepars, opt_param
%-RasterCell    Raster trials number of entries
%  
%OUTPUT: RASTER_LOGPROB
% .perBin       vector ofthe logprob of each bin, collapse out trials
%               each trial of Bin independent event..simple summation
% .perTrial     cell with trials entries
%               each Bin of each trial ind event.. simple summation
% ._avgProbability
%               additive avg of LogProbs then exponentiated
%               good measure of avg probability of each individiual bin
function [Raster_LogProb]     = GLMcompute_rasterLogProbwrapper(Basepars,SimandMetrics, p_opt)

GLMPars                       = Basepars.GLMPars;
[StimulusPars DirPars ]       = Directories_Params_func_all(Basepars.exp_nm, GLMPars.debug);
MetricPars = SimandMetrics.MetricPars; count = 1;
%%%%%%%%%% Directory Fixing %%%%%%%%%%%
fit_type = GLMPars.fit_type; Head_CTYPE = Basepars.celltype;

%%
for count = 1 :2
    %%
	if count == 1,  stim_vrf = 'BW'  ;  Slv_StimPars = StimulusPars.BW  ;      end
    if count == 2,  stim_vrf = 'NSEM';  Slv_StimPars = StimulusPars.NSEM;      end
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
        movie_width  = StimulusPars.NSEM.Width;
        movie_height = StimulusPars.NSEM.Height;
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
        [lcif_cp, lcif_perneighbor] = logintensityfunction_coupling (p_opt, Basepars, stim_vrf, GLMPars.debug); 
    end
    
    CIFcomponents.lcif_mu = lcif_mu;
    CIFcomponents.lcif_kx = lcif_kx;
    if  Basepars.Coupling  && ~isempty(Basepars.cp_Neighbors) 
        CIFcomponents.lcif_cp = lcif_cp;
    end
    CIFcomponents.cif_psgain = ps_cifgain;   
    clear lcif_cp lcif_mu lcif_kx p_ps psfilter ps_cifgain logical_BiDirect lcif_perneighbor
    display('ModelLoaded into CIFcomponents structure')
    
    %% Choose Verification Dateset and load Rasters
    Cell_ID = Basepars.headid; expnm_string = Basepars.exp_nm; logical_debug = GLMPars.debug; rastertype_string = stim_vrf;
    bins_perstimframe =  Basepars.spikebins_perstimframe;  % large enough so atmost  one spike per bin
    plot_logical = false; 
    [~, Raster_sptimeSecs ] = load_rasters(Cell_ID, expnm_string, rastertype_string, bins_perstimframe, logical_debug, plot_logical);
    clear bins_perstimframe debug rastertype_string expnm_string Cell_ID  logical_debug 
    
    %%
	if Basepars.Coupling  && ~isempty(Basepars.cp_Neighbors)
        blocks = max(size(CIFcomponents.lcif_cp));
    else 
        blocks = max(size(Slv_StimPars.Raster_Blocks));
	end
    clear RasterCell
    RasterCell = cell(blocks, 1);
    for i_block = 1 : blocks
        RasterCell{i_block}.block  = Raster_sptimeSecs{i_block}.block;
        RasterCell{i_block}.spikes = Raster_sptimeSecs{i_block}.spikes;
    end
    if count == 1, duration = Slv_StimPars.nsec_o;   end
    if count == 2, duration = Slv_StimPars.nsec_o;   end

    Raster_LogProb = compute_rasterLogProb(RasterCell , CIFcomponents, duration);
    
    if count == 1, BW_LogProb   = Raster_LogProb; clear Raster_LogProb; end
    if count == 2, NSEM_LogProb = Raster_LogProb; clear Raster_LogProb; end
end
Raster_LogProb.BW   =   BW_LogProb;
Raster_LogProb.NSEM = NSEM_LogProb; 