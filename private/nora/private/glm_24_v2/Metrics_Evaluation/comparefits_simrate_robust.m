% AkHeitman 2014-03-8   generate a solid metric
% Fix on a fit_type comparison that we trust
% Deal with Simulating the GLM later .. 2-spikes?
% Later incorporate the lag in test seconds

clear; close all
GLMdir ='/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/GLM' ;
basedir = '/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/GLM/Metric_Testing/';
        map_type = 'mapPRJ';
filtertype = 'fixedSP';
%exp_nm = '2012-09-27-3'; expname = 'expB'; fit_type = 'NSEM';  %cmodeltype
%= '8pix_Identity_8pix';
shortlist = true;
binnedrast_dur0 = '/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/GLM/Rasters/binnedspikes_logprob_bpdf10/';
cmodel = '8pix_Identity_8pix';  
version = 'BWBW_vs_BWNSEM'

i_exp = 1; i_test = 2;

%eval(sprintf('load %s/%s_logprobs_shortlist.mat' , binnedrast_dur0, test_type));



%i_test = 1;  i_exp = 1;
%savefile= 
%%
for i_exp = 1:2
    %%
    if shortlist
        if i_exp == 1, exp_nm = '2012-08-09-3' ; expname = 'expA'; head_cellID = [1471  3676 5086 5161 841 1426 1772 2101 1276];
        if i_exp == 2, exp_nm = '2012-09-27-3';  expname = 'expB'; head_cellID = [1 31 301 1201 1726 91 1909 2360 6858];  end
        if i_exp == 3, exp_nm = '2013-08-19-6';  expname = 'expC'; head_cellID = [ 2824 3167 3996 5660 6799]; end % OFFPar doesn't converge for NSEM [737 1328 1341 2959  5447];
        if i_exp == 4, exp_nm = '2013-10-10-0';  expname = 'expD'; head_cellID = [32 768 2778 4354 5866 7036];end % OFFPar doesn't converge for NSEM [346 1233 3137  5042 5418];
    end

    compareWNNSEM.cell_ids  = head_cellID;   
    compareWNNSEM.note      = 'Row 1 is WN-WN Row 2 is NSEM-NSEM';
    compareWNNSEM.note2     = version;
    compareWNNSEM.note3     = 'For comparing fits ,  NSEM simulations do not always converge,   ONPar 2012-08-09-3 and OFFPar fpr 2013-08-19-6 2013-10-10-0'; 
    compareWNNSEM.ratediff_normrate = zeros(2, length(head_cellID));
    compareWNNSEM.ratediff_normvar = zeros(2, length(head_cellID));
    compareWNNSEM.ratediff_normrate5 = zeros(2, length(head_cellID));
    compareWNNSEM.ratediff_normvar5 = zeros(2, length(head_cellID));    
    compareWNNSEM.ratediff_normrate10 = zeros(2, length(head_cellID));
    compareWNNSEM.ratediff_normvar10 = zeros(2, length(head_cellID));
    compareWNNSEM.avgfiringrate   =  zeros(2,length(head_cellID) );
    compareWNNSEM.ONPar           =  zeros(1,length(head_cellID) ); 
    compareWNNSEM.OFFPar          =  zeros(1,length(head_cellID) ); 
    
    compareWNNSEM.logprob_glm_raw = zeros(2,length(head_cellID) );
    compareWNNSEM.logprob_mr_rwa  = zeros(2,length(head_cellID) );
    compareWNNSEM.logprob_glm_bps = zeros(2,length(head_cellID) );
    compareWNNSEM.logprob_mr_bps  = zeros(2,length(head_cellID) );
    
      
   %% 
    for i_test = 1:2
        %%
        if i_test  == 1
            if strcmp(version , 'BWBW_vs_BWNSEM')
                fit_type = 'BW';
                test_type = 'BW';
            end    
        end
        if i_test  == 2
            if strcmp(version , 'BWBW_vs_BWNSEM')
                fit_type = 'BW';
                test_type = 'NSEM';
            end    
        end
        [StimulusPars DirPars datarun_slv datarun_mas] = Directories_Params_v19_split(exp_nm, fit_type, 'mapPRJ');
        if strcmp(test_type , 'NSEM')
            SPars = StimulusPars.NSEM;
        elseif strcmp(test_type, 'BW')
            SPars = StimulusPars.BW;
        end
        clear DirPars StimulusPars
        
        %%
        binnedrast_dur = sprintf('%s/%s_%s', binnedrast_dur0 , exp_nm , test_type);         
        orgspikesdir = sprintf('/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/GLM/%s/%s_mapPRJ/BlockSpikesandRasters', exp_nm, test_type);  
        stimulidir = '/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli';
        [blockedmoviecell,  novelmoviestats, origmatfile] = loadmoviematfile(exp_nm , test_type, cmodel, 'testmovie')   
        testmovie = blockedmoviecell{1};

        fittype = sprintf('%s_%s/%s_Fit_preCosyne' , fit_type, map_type, fit_type);
        fitparams  = sprintf('%s_ps20_cpOFF/bin10_blk55_tolfun6/%s', filtertype, cmodel);
        if strcmp(exp_nm, '2013-10-10-0')
            fitparams =  sprintf('%s_ps20_cpOFF/bin10_blk27_tolfun6/%s', filtertype, cmodel);
        end
        if strcmp(fit_type, 'BW')
             fitparams = sprintf('%s_ps20_cpOFF/bin10_blk50_tolfun6/%s', filtertype, cmodel);
        end
        if strcmp(fit_type, 'BW') && strcmp(filtertype , 'rk2')
             fitparams = sprintf('%s_ps20_cpOFF/bin10_blk50_tolfun3/%s', filtertype, cmodel);
        end
        savecomparisonsdir = sprintf('%s/%s/%s_ps20_cpOFF', basedir , cmodel,filtertype); 
        modeltestdir  = sprintf('%s/%s/%s_ps20_cpOFF/%s_test/%s_fit/%s', basedir , cmodel,filtertype,test_type,fit_type,exp_nm);
        if ~exist(modeltestdir), mkdir(modeltestdir); end
        cid = head_cellID(1)
%%

for cid = head_cellID
    %%  Load up Basepars,  binned spies etc.  
    index = find( head_cellID == cid);
        if ~isempty(find(datarun_mas.cell_types{1}.cell_ids == cid))
            celltype = 'ONPar';
             compareWNNSEM.ONPar(index) = 1;
            
        end
        if ~isempty(find(datarun_mas.cell_types{2}.cell_ids == cid))
            celltype = 'OFFPar';
            compareWNNSEM.OFFPar(index) = 1;
        end

    sname = sprintf('%s_%d', celltype, cid)
    if strcmp(filtertype, 'rk2') && strcmp(fit_type, 'NSEM')
        load(sprintf('%s/%s/%s/%s/L13_%s.mat' , GLMdir, exp_nm,fittype, fitparams , sname));
    end
    
    %% comebackheresoon!! %% 
        
    load(sprintf('%s/%s/%s/%s/%s.mat' , GLMdir, exp_nm,fittype, fitparams , sname));
    eval(sprintf('load %s/spikerates_binnedlogprobs_%s.mat', binnedrast_dur , sname))
    eval(sprintf('load %s/organizedspikes_%s.mat', orgspikesdir , sname))
%%
    testmovie_ROI = testmovie.matrix(Basepars.ROI.xdim , Basepars.ROI.ydim, :);
    testmovie_ROI = double(testmovie_ROI);
    testmovie_ROI = testmovie_ROI / max(testmovie_ROI(:));
    testmovie_ROI = testmovie_ROI - inputstats.mu_avgIperpix;

    
    display(sprintf('~~~ Minimum stim is %d ---' , min(testmovie_ROI(:))));
    display(sprintf('~~~ Maximum stim is %d ---' , max(testmovie_ROI(:))));
    

    klen = length(Basepars.ROI.xdim); frames = SPars.raster_frames;
    testmovie_ROI = reshape(testmovie_ROI , [klen^2,frames]);
    figure; hist(testmovie_ROI(:), 20); title('testmovie distribution')

    %% Set up CIF Components
    p_opt = Basepars.p_opt; ps_basis = Basepars.ps_basis;
    dt        = Basepars.dt; tstim =Basepars.tstim;
    parInd    = Basepars.paramind;
    MU = p_opt(parInd.MU);
    PS = ps_basis * p_opt(parInd.PS);

    if strcmp(Basepars.k_filtermode , 'fixedSP')
        spfilter = Basepars.spfilter;
        if size(spfilter,2) > size(spfilter,1), spfilter = spfilter'; end
        K  = (spfilter)*(p_opt(parInd.L)');
        timefilter = p_opt(parInd.L)';
        testmovie_scalar = (spfilter') * testmovie_ROI;
        lcif_kx_frame = fastconvAH( testmovie_scalar, timefilter , 1, frames,0);     
        %figure; subplot(2,1,1); plot( timefilter); subplot(2,1,2); plot(lcif_kx_frame);
    end
    
    if strcmp(Basepars.k_filtermode , 'rk2')      
        Z = Basepars.paramind;
        p_opt = Basepars.opt_param.p;
        MU = p_opt(Z.MU);
        K1 = p_opt(Z.SPACE1)*(p_opt(Z.TIME1)');
        K2 = p_opt(Z.SPACE2)*(p_opt(Z.TIME2)');
        K = K1 + K2;        
        lcif_kx_frame = zeros(size(testmovie_ROI));
        for i_spot = 1:klen^2
            lcif_kx_frame(i_spot,:) = fastconvAH( testmovie_ROI(i_spot,:), K(i_spot,:) , 1, frames,0);
        end
        lcif_kx_frame = sum(lcif_kx_frame);
    end

    display('binning the lcif components')

    sbpf   = Basepars.spikebins_perstimframe;
    lcif_kx = zeros(frames*sbpf , 1);
    for i = 1 : frames
            indices = 1 + (i-1)*sbpf: i*sbpf ;
            lcif_kx(indices) = lcif_kx_frame(i);
    end
    lcif_mu = MU * ones (frames * sbpf,1); 

    CIFcomponents.lcif_mu = lcif_mu;
    CIFcomponents.lcif_kx = lcif_kx;
    CIFcomponents.lcif_kxframe        = lcif_kx_frame; % unnnecesarry arg for simulateandmodel
   % CIFcomponents.scalar_stimframe    =  testmovie_scalar; % unnnecesarry arg for simulateandmodel
    CIFcomponents.cif_psgain = exp(PS);


    RasterBlocks = SPars.FitBlocks - 1;
    SimPars.seconds   = SPars.nsec_o;
    SimPars.bin_dt    = dt;
    SimPars.PostSpike = true;
    SimPars.Coupling = false;
    SimPars.trials   = 1;
    SimPars.simulations_pertrial = length(RasterBlocks);
    %Simulations = [];
    %
    Simulations = simulate_model(SimPars , CIFcomponents);
    %



    
    %% Simulation Rates
     binnedspRast = double(spikerates_binnedlogprobs.binnedspikes);
     bins = size(binnedspRast,2); 
    
     blks = length(RasterBlocks);
     binnedspGLM = zeros(blks, bins);
     for i_blk = 1:blks
         binnedspGLM(i_blk , : ) = double( Simulations{i_blk}.binnedlogical );
     end
     
     rateGLM = sum(binnedspGLM ,1 ) / blks;     
     rateRast      = sum(binnedspRast ,1 ) / blks;
     
     rateGLM5 = sum(reshape(rateGLM, [5 , bins/5]));
     rateRast5 = sum(reshape(rateRast, [5 , bins/5]));
     
     
     rateGLM10 = sum(reshape(rateGLM, [10 , bins/10]));
     rateRast10 = sum(reshape(rateRast, [10 , bins/10]));
     
     rateGLM20 = sum(reshape(rateGLM, [20 , bins/20]));
     rateRast20 = sum(reshape(rateRast, [20 , bins/20]));
     
     
     c_rate_normrastrate  =   sum((rateGLM - rateRast).^2) / sum((rateRast).^2)
     c_rate_normrastvar   =   sum((rateGLM - rateRast).^2) / sum( (rateRast-mean(rateRast)).^2 )
     c_rate_normrastrate5 =   sum((rateGLM5 - rateRast5).^2) / sum((rateRast5).^2)
     c_rate_normrastvar5  =   sum((rateGLM5 - rateRast5).^2) / sum( (rateRast5-mean(rateRast5)).^2 )
     c_rate_normrastrate10 =   sum((rateGLM10 - rateRast10).^2) / sum((rateRast10).^2)
     c_rate_normrastvar10  =   sum((rateGLM10 - rateRast10).^2) / sum( (rateRast10-mean(rateRast10)).^2 )
     c_rate_normrastrate20 =   sum((rateGLM20 - rateRast20).^2) / sum((rateRast20).^2)
     c_rate_normrastvar20  =   sum((rateGLM20 - rateRast20).^2) / sum( (rateRast20-mean(rateRast20)).^2 )
     
     
    
     
    
     compareWNNSEM.ratediff_normrate(i_test, index)   = c_rate_normrastrate;
     compareWNNSEM.ratediff_normvar(i_test, index)    = c_rate_normrastvar ;
     compareWNNSEM.ratediff_normrate5(i_test, index)  = c_rate_normrastrate5;
     compareWNNSEM.ratediff_normvar5(i_test, index)   = c_rate_normrastvar5; 
     compareWNNSEM.ratediff_normrate10(i_test, index) = c_rate_normrastrate10;
     compareWNNSEM.ratediff_normvar10(i_test, index)  = c_rate_normrastvar10  
     compareWNNSEM.ratediff_normrate20(i_test, index) = c_rate_normrastrate20;
     compareWNNSEM.ratediff_normvar20(i_test, index)  = c_rate_normrastvar20  ;     
     compareWNNSEM.avgfiringrate(i_test, index)       = sum(binnedspRast(:)) / (blks * SPars.nsec_o) ;

  
     
     
    %% Compute Log Prob
   %{
    CIF = CIFcomponents
    CIF.lcif_ps = log(CIF.cif_psgain);
    lcif_ps = log(CIF.cif_psgain);
    base_lcif = (CIF.lcif_mu + CIF.lcif_kx)';
    lastbin = bins - length(CIF.lcif_ps)
    rast_lcif0 = repmat(base_lcif , [blks , 1]);
    convolvedps = fastconvAH(binnedspRast , [0 lcif_ps'],blks, bins);
  
    rast_lcif = rast_lcif0 + convolvedps;
    rast_cif = exp(rast_lcif);
    prob_spike   = Basepars.dt * rast_cif;
    toohighind = find(prob_spike >= 1);
    if ~isempty(toohighind)
        prob_spike(toohighind) = 1 - eps;
    end
    prob_nospike = 1 - prob_spike;
    rawbinprob = prob_nospike;
    for i_blk = 1 : length(RasterBlocks)
        spikebin = find(binsp( i_blk, :) );
        rawbinprob(i_blk , spikebin ) = prob_spike(spikebin) ;
    end    
    logbinprob = log(rawbinprob);
    
    
    %{   
    lagsecs   = 1;
    lagframes = 120*lagsecs;
    testframes = SPars.raster_frames - lagframes;
    %}    
    logtrialprob        = sum( logbinprob , 2 );
    avglogprob          = mean(logtrialprob);
    avglogprobperframe  = avglogprob / testframes;
    avgframeprob        = exp(avglogprobperframe)
    normalizedframeprob = avgframeprob / exp(spikerates_binnedlogprobs.perframelogprob)


    spratenullogprob = ((SPars.nsec_o-lagsecs) / SPars.nsec_o )   * spikerates_binnedlogprobs.spratenullogprob;
    onepercentdur    =        spratenullogprob / log(.01);        
    betterthanchance_onepercent = ( exp( avglogprob /  onepercentdur) ) / .01 
    normalized_btc = betterthanchance_onepercent / spikerates_binnedlogprobs.betterthanchance_onepercent



    compareWNNSEM.logprobs(index)          =  avglogprob ; 
    compareWNNSEM.perframeprobs(index)     =  avgframeprob;
    compareWNNSEM.nperframeprobs(index)    =  normalizedframeprob;
    compareWNNSEM.avgfiringrate(index)     =  mean(spikerates_binnedlogprobs.framerate_hz);
    compareWNNSEM.btc(index)               =  betterthanchance_onepercent; 
    compareWNNSEM.nbtc(index)              =  normalized_btc; 

%}




end

    end
    eval(sprintf('save %s/compareperformance_%s_%s.mat compareWNNSEM', savecomparisonsdir, version,exp_nm));
end

    

