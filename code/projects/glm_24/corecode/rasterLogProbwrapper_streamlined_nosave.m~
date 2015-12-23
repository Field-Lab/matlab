
% Dirty Code for EJ's figure  2013-02-22
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
function [Raster_LogProb] = rasterLogProbwrapper_streamlined_nosave(Basepars, SimandMetrics ,p_opt)
 
[StimulusPars DirPars datarun_all] = Directories_Params_func_all(Basepars.exp_nm, false)
%DirPars = Basepars.DirPars;
GLMPars              = Basepars.GLMPars;
datarun_BW{1}  = datarun_all.master;  datarun_BW{2}   = datarun_all.BW_slv;
datarun_NSEM{1} = datarun_all.master; datarun_NSEM{2} = datarun_all.NSEM_slv; 


GLMPars                       = Basepars.GLMPars;

%MetricPars = SimandMetrics.MetricPars; count = 1;
%%%%%%%%%% Directory Fixing %%%%%%%%%%%
fit_type = GLMPars.fit_type; Head_CTYPE = Basepars.celltype;

%%
for count = 1 :2
    %%
	if count == 1,  stim_vrf = 'BW'  ;  Slv_StimPars = StimulusPars.BW  ;   datarun = datarun_BW;     end
    if count == 2,  stim_vrf = 'NSEM';  Slv_StimPars = StimulusPars.NSEM;   datarun = datarun_NSEM;   end
	%%  LOAD THE MOVIES AND TESTRUN 
    
    %%%%  LOADING TESTRUN %%%%%%%%%
    load(sprintf('%s/%s/%s/prep_STAandROI/testrun_id%d_%dsq',DirPars.output_dir,'BW',Basepars.celltype,Basepars.headid,GLMPars.K_slen));  
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
    
    
    if isfield(Basepars.paramind, 'FR')   
         mov       = testmovie_ROI;
        t_filter = ones(size(Basepars.STA,1),1) * Basepars.int_time_filter'; 

        %%% m_coit   = movie convolved over the intergration time
        m_coit =  fastconv(mov,t_filter,Basepars.k_spacepixels,size(mov,2),Basepars.padval);
        pos_parts = find(m_coit >= 0) ;
        neg_parts = find(m_coit < 0 ) ;
        if strcmp(Basepars.celltype, 'ON-Parasol'), rec_m_coit = m_coit; rec_m_coit(neg_parts) = 0; end     
        if strcmp(Basepars.celltype, 'OFF-Parasol'), rec_m_coit = m_coit; rec_m_coit(pos_parts) = 0; end
        
        sqrec = rec_m_coit .^2;
        
        if strcmp(Basepars.celltype, 'OFF-Parasol'), sqrec = - sqrec; end
        
        
        
       % sqrec = Basepars.rect_full ;
        lcif_rec0 =(sqrec') * p_opt(Basepars.paramind.FR);
        lcif_rec  = reprows(lcif_rec0, Basepars.spikebins_perstimframe);
        CIFcomponents.lcif_rec = lcif_rec;
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