% AKHeitman 2014-03-23
% Effectively Done as of 2014-03-23
% Making Robust Model Comparison 
% Decided to fix on  (LogProb model - LogProb null) / spikes
% Normalized by (LogProb UOP - LogProb null) / spikes



%%%% Set up directories %%%%
clear; close all; clear all
GLMdir ='/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/GLM' ;
modelcompare_basedir = '/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/GLM/Model_Compare/';
shortlist = true;

%%% Modulation and mod_version %%
fixed.fitversion = 'Fit';
% define hatges changed and what stays fixed
% modulation is either, "fit_and_test" or "model_parameters" or "other"
% "fit_and_test" referes to NSEM or WN for fitting and testing
% "model_parameters"  should be flexible for any model change

%{
modulation  = 'model_parameters';
mod_version = 'NSEM_fixedSPvsrk2_8pix_ID_8pix'; purpose = 'Define a benchmark NSEM improvement level: the switch to rk2 filters from fixed spatial filters.  No cone model';
comparenames{1} = 'fixedsp';
comparenames{2} = 'rk2';
comp_params{1}.filtertype = 'fixedSP';
comp_params{2}.filtertype = 'rk2';
comp_params{1}.cmodel    ='8pix_Identity_8pix';
comp_params{2}.cmodel    ='8pix_Identity_8pix';
fixed.fit_type = 'NSEM'; fixed.test_type = 'NSEM';
%}

%
modulation  = 'model_parameters';
mod_version = 'WN_fixedSPvsrk2_8pix_ID_8pix'; purpose = 'Define a benchmark WN improvement level: the switch to rk2 filters from fixed spatial filters. No cone model';
comparenames{1} = 'fixedsp';
comparenames{2} = 'rk2';
comp_params{1}.filtertype = 'fixedSP';
comp_params{2}.filtertype = 'rk2';
comp_params{1}.cmodel    ='8pix_Identity_8pix';
comp_params{2}.cmodel    ='8pix_Identity_8pix';
fixed.fit_type = 'BW'; fixed.test_type = 'BW';
%}



%{
modulation  = 'model_parameters';
mod_version = 'WN_ONFFvsLin_8pix_1e4_8pix'; purpose = 'Examine effect of ONOff Hard rectification channels in WN';
comparenames{1} = 'single linear channel';
comparenames{2} = 'hard rect ON and Off channel';
comp_params{1}.filtertype = 'fixedSP';
comp_params{2}.filtertype = 'OnOff_hardrect_fixedSP_STA';
comp_params{1}.cmodel    ='8pix_Model1_1e4_8pix';
comp_params{2}.cmodel    ='8pix_Model1_1e4_8pix';
fixed.fit_type = 'BW'; fixed.test_type = 'BW';
%}
%{
modulation  = 'model_parameters';
mod_version = 'NSEM_ONFFvsLin_8pix_1e4_8pix'; purpose = 'Examine effect of ONOff Hard rectification channels in NSEM';
comparenames{1} = 'single linear channel';
comparenames{2} = 'hard rect ON and Off channel';
comp_params{1}.filtertype = 'fixedSP';
comp_params{2}.filtertype = 'OnOff_hardrect_fixedSP_STA';
comp_params{1}.cmodel    ='8pix_Model1_1e4_8pix';
comp_params{2}.cmodel    ='8pix_Model1_1e4_8pix';
fixed.fit_type = 'NSEM'; fixed.test_type = 'NSEM';
%}


%{
modulation  = 'model_parameters';
mod_version = 'NSEM_ONFFvsLin_8pix_ID_8pix'; purpose = 'Examine effect of ONOff Hard rectification channels in NSEM, no cone model';
comparenames{1} = 'single linear channel';
comparenames{2} = 'hard rect ON and Off channel';
comp_params{1}.filtertype = 'fixedSP';
comp_params{2}.filtertype = 'OnOff_hardrect_fixedSP_STA';
comp_params{1}.cmodel    ='8pix_Identity_8pix';
comp_params{2}.cmodel    ='8pix_Identity_8pix';
fixed.fit_type = 'NSEM'; fixed.test_type = 'NSEM';
%}
%}
%%% Other modulation options
%{





%modulation = 'fit_and_test';
mod_version = 'WNWN_vs_NSEMNSEM'; purpose = 'Show WN GLM fits better than NSEM GLM fits'; 
comparenames{1} = 'fit/test by WN'; 
comparenames{2} = 'fit/test by NSEM';

%mod_version = 'WNWN_vs_WNNSEM'; purpose = 'Show WN fits do not generalize to NSEM'; 
%comparenames{1} = 'fit/test by WN'; 
%comparenames{2} = 'fit with WN, test with NSEM';

%mod_version = 'NSEMWN_vs_NSEMNSEM'; purpose = 'Show NSEM fits do not generalize to NSEM'; 
%comparenames{1} = 'fit with NSEM, test with WN';
%comparenames{2} = 'fit/test by NSEM'; 
modulation  = 'model_parameters';
mod_version = 'WN_conemodel_8pix_1e4_8pix'; purpose = 'Show Cone model has minimal effect for WN fits';
comparenames{1} = 'fit without cone model';
comparenames{2} = 'fit with cone model 8p1e48p';
comp_params{1}.filtertype = 'fixedSP';
comp_params{2}.filtertype = 'fixedSP';
comp_params{1}.cmodel    ='8pix_Identity_8pix';
comp_params{2}.cmodel    ='8pix_Model1_1e4_8pix';
fixed.fit_type = 'BW'; fixed.test_type = 'BW';
%}

fixed.map_type = 'mapPRJ';
if ~strcmp(modulation, 'model_parameters')
    fixed.filtertype ='fixedSP';
    fixed.cmodel = '8pix_Identity_8pix';
end


%%% Write up savedir creat saving structure %%%
savecomparisonsdir = sprintf('%s/%s/%s', modelcompare_basedir , modulation, mod_version);
if ~exist(savecomparisonsdir, 'dir'), mkdir(savecomparisonsdir); end
full = cell(4,1);
i_exp = 1; i_test = 1;
%%
for i_exp = 1:4
    %% Setup model_comparison structure  one per dataset
    if shortlist
        if i_exp == 1, exp_nm = '2012-08-09-3' ; expname = 'expA'; head_cellID = [1471  3676 5086 5161 841 1426 1772 2101 1276]; end
        if i_exp == 2, exp_nm = '2012-09-27-3';  expname = 'expB'; head_cellID = [1 31 301 1201 1726 91 1909 2360 6858];  end
        if i_exp == 3, exp_nm = '2013-08-19-6';  expname = 'expC'; head_cellID = [ 2824 3167 3996 5660 6799 737 1328 1341 2959  5447]; end % OFFPar doesn't converge for NSEM [737 1328 1341 2959  5447];
        if i_exp == 4, exp_nm = '2013-10-10-0';  expname = 'expD'; head_cellID = [32 768 2778 4354 5866 7036 346 1233 3137  5042 5418];end % OFFPar doesn't converge for NSEM [346 1233 3137  5042 5418];
    end
    clear model_comparison
    model_comparison.exp_nm            = exp_nm;
    model_comparison.cell_ids          = head_cellID;   
    model_comparison.comparison_type   = modulation;
    model_comparison.mod_version       = mod_version;
    model_comparison.purpose           = purpose;
    model_comparison.comparenames      = comparenames;
    model_comparison.savedir           = savecomparisonsdir;
    model_comparison.note0             = 'Logarithmic Probability based comparisons'; 
    model_comparison.note1             = 'UOP: Uncoditioned Optimal Perfomance, mean rate model'; 
    model_comparison.note2             = 'BPS: bits per spike'; 
    model_comparison.note3             = 'Bits stand for LogProb - LogProb of uniform rate model';
    model_comparison.note4             = 'Uniform rate model represents no information but total spikes';
    model_comparison.ONPar                = zeros(1,length(head_cellID)); 
    model_comparison.OFFPar               = zeros(1,length(head_cellID)); 
    model_comparison.logprob_glm_raw      = zeros(2,length(head_cellID));
    model_comparison.logprob_glm_bpspike  = zeros(2,length(head_cellID));
    model_comparison.logprob_glm_bpsec    = zeros(2,length(head_cellID));
    model_comparison.logprob_uop_raw      = zeros(2,length(head_cellID));
    model_comparison.logprob_uop_bpspike  = zeros(2,length(head_cellID));
    model_comparison.logprob_uop_bpsec    = zeros(2,length(head_cellID));
    
    for i_test = 1:2
        %% Fit / test type and model params
        if strcmp(modulation, 'fit_and_test')
            if i_test  == 1
                if strcmp(mod_version , 'WNWN_vs_WNNSEM'),      fit_type = 'BW';    test_type = 'BW'; end
                if strcmp(mod_version , 'WNWN_vs_NSEMNSEM'),    fit_type = 'BW';    test_type = 'BW'; end
                if strcmp(mod_version , 'NSEMWN_vs_NSEMNSEM'),  fit_type = 'NSEM';  test_type = 'BW'; end
            end
            if i_test  == 2
                if strcmp(mod_version , 'WNWN_vs_WNNSEM'),     fit_type = 'BW';     test_type = 'NSEM'; end
                if strcmp(mod_version , 'WNWN_vs_NSEMNSEM'),   fit_type = 'NSEM';   test_type = 'NSEM'; end
                if strcmp(mod_version , 'NSEMWN_vs_NSEMNSEM'), fit_type = 'NSEM';   test_type = 'NSEM'; end
            end
        else
            fit_type = fixed.fit_type; test_type = fixed.test_type;
        end        
        if strcmp(modulation, 'model_parameters')
            filtertype = comp_params{i_test}.filtertype;
            cmodel     = comp_params{i_test}.cmodel;
        else
            cmodel = fixed.cmodel; 
            filtertype = fixed.filtertype;
        end
        map_type = fixed.map_type    
        fittype = sprintf('%s_%s/%s_%s' , fit_type, map_type, fit_type, fixed.fitversion);  

        [StimulusPars DirPars datarun_slv datarun_mas] = Directories_Params_v20(exp_nm, fit_type, 'mapPRJ');
        if strcmp(test_type , 'NSEM')
            SPars = StimulusPars.NSEM;
        elseif strcmp(test_type, 'BW')
            SPars = StimulusPars.BW;
        end
        orgspikesdir = sprintf('/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/GLM/%s/%s_mapPRJ/BlockSpikesandRasters', exp_nm, test_type);  
        stimulidir = '/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli';
        [blockedmoviecell,  novelmoviestats, origmatfile] = loadmoviematfile(exp_nm , test_type, cmodel, 'testmovie')   
        testmovie = blockedmoviecell{1};
        clear DirPars StimulusPars stimulidir origmatfile blockedmoviecell              
        
        if strcmp(fit_type, 'BW')
             fitparams = sprintf('%s_ps20_cpOFF/bin10_blk50_tolfun5/%s', filtertype, cmodel);
        end
        if strcmp(fit_type, 'NSEM')
             fitparams = sprintf('%s_ps20_cpOFF/bin10_blk55_tolfun5/%s', filtertype, cmodel);
        end
        if strcmp(fit_type, 'BW') && strcmp(filtertype , 'rk2')
             fitparams = sprintf('%s_ps20_cpOFF/bin10_blk50_tolfun5/%s', filtertype, cmodel);
        end
        if strcmp(fit_type, 'NSEM') && strcmp(filtertype , 'rk2')
             fitparams = sprintf('%s_ps20_cpOFF/bin10_blk55_tolfun5/%s', filtertype, cmodel);
        end
        if strcmp(exp_nm, '2013-10-10-0') && strcmp( fit_type , 'NSEM')
            fitparams =  sprintf('%s_ps20_cpOFF/bin10_blk27_tolfun5/%s', filtertype, cmodel);
        end
        cid = head_cellID(1)

%% Individual Cell Loops . try to make this a functional unit for evaluation w/n GLM fitting
for cid = head_cellID
    %%%% Load up Basepars,  binned spikes , cell specific information  %%%
    index = find( head_cellID == cid);
	if ~isempty(find(datarun_mas.cell_types{1}.cell_ids == cid))
            celltype = 'ONPar';
             model_comparison.ONPar(index) = 1;
            
    end
	if ~isempty(find(datarun_mas.cell_types{2}.cell_ids == cid))
            celltype = 'OFFPar';
            model_comparison.OFFPar(index) = 1;
    end
    params.sname = sprintf('%s_%d', celltype, cid)
    load(sprintf('%s/%s/%s/%s/%s.mat' , GLMdir, exp_nm,fittype, fitparams , params.sname));
    eval(sprintf('load %s/organizedspikes_%s.mat', orgspikesdir , params.sname))
    
    
    
    params.tstim = Basepars.tstim;
    params.dt    = Basepars.dt;
    params.bindur = params.tstim / Basepars.spikebins_perstimframe; 
    params.bins = Basepars.spikebins_perstimframe *length(SPars.testframes);  
    params.evalblocks = SPars.evalmodel_Blocks;
	params.trials = length(params.evalblocks);  
    params.frames = length(SPars.testframes);
    params.klen = length(Basepars.ROI.xdim);
    params.testdur_seconds = max(SPars.testseconds) - min(SPars.testseconds) ;  
    
	logicalspike = zeros( length(params.evalblocks) , params.bins) ;         
	for i_blk = 1 : length(params.evalblocks)
            blknum = params.evalblocks(i_blk);
            sptimes = organizedspikes.block.t_sp_withinblock{blknum} - SPars.testseconds(1);
            sptimes = sptimes(find(sptimes > 0 ) );
            for i_sp = 1:length(sptimes)
                spt = sptimes(i_sp);
                binnumber = ceil(spt / params.bindur );
                logicalspike( i_blk, binnumber )  =  logicalspike( i_blk,binnumber ) + 1;
            end
    end    
    clear tstim bindur evalblocks testbins sptimes spt binnumber 
    
    testmovie_ROI = testmovie.matrix(Basepars.ROI.xdim , Basepars.ROI.ydim, :);
    testmovie_ROI = double(testmovie_ROI);
    testmovie_ROI = testmovie_ROI / max(testmovie_ROI(:)); 
    testmovie_ROI = testmovie_ROI - novelmoviestats.mu_avgIperpix;
    testmovie_ROI = testmovie_ROI(:,:, SPars.testframes);
    testmovie_ROI = reshape(testmovie_ROI , [params.klen^2,params.frames]);
    
    display(sprintf('~~~ Minimum stim is %d ---' , min(testmovie_ROI(:))));
    display(sprintf('~~~ Maximum stim is %d ---' , max(testmovie_ROI(:))));   
    clf; hist(testmovie_ROI(:), 20); title('testmovie distribution')
    %% Set up CIF Components
    clear lcif_kx_frame
    
    p_opt     = Basepars.p_opt; ps_basis = Basepars.ps_basis;
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
        lcif_kx_frame = fastconvAH( testmovie_scalar, timefilter , 1, length(testmovie_scalar),0);     
        %figure; subplot(2,1,1); plot( timefilter); subplot(2,1,2); plot(lcif_kx_frame);
    end
    
    if strcmp(Basepars.k_filtermode , 'OnOff_hardrect_fixedSP_STA')
        spfilter_pos = Basepars.spfilter_pos;
        spfilter_neg = Basepars.spfilter_neg;

        if size(spfilter_pos,2) > size(spfilter_pos,1), spfilter_pos = spfilter_pos';  spfilter_neg = spfilter_neg';end
    %	K1 = (spfilter)*(p_opt(parInd.TIME1)');
        timefilter1 = p_opt(parInd.TIME1)';
        timefilter2 = p_opt(parInd.TIME2)';

        testmovie_ROIpos = testmovie_ROI;
        testmovie_ROIneg = testmovie_ROI;
        negind = find(testmovie_ROI(:) < 0 );
        posind = find(testmovie_ROI(:) > 0 );

        testmovie_ROIpos(negind) = 0;
        testmovie_ROIneg(posind) = 0;

        testmovie_pos = (spfilter_pos') * testmovie_ROIpos;
        testmovie_neg = (spfilter_neg') * testmovie_ROIneg;

        lcif_kx_frame1 = fastconvAH( testmovie_pos, timefilter1 , 1, length(testmovie_pos),0);
        lcif_kx_frame2 = fastconvAH( testmovie_neg, timefilter2 , 1, length(testmovie_neg),0); 
        lcif_kx_frame  = lcif_kx_frame1 + lcif_kx_frame2;
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
        for i_spot = 1:Basepars.k_spacepixels
            lcif_kx_frame(i_spot,:) = fastconvAH( testmovie_ROI(i_spot,:), K(i_spot,:) , 1, length(testmovie_scalar),0);
        end
        lcif_kx_frame = sum(lcif_kx_frame);
    end

    display('binning the lcif components')
    
    lcif_kx0 = reshape( repmat(lcif_kx_frame, 10, 1) , 1 , params.bins);
    lcif_mu0 = MU * ones (1,params.bins); 
    
    lcif_mu = repmat(lcif_mu0 , params.trials, 1);
    lcif_kx = repmat(lcif_kx0 , params.trials, 1);
    
    clear sbpf;   
    lcif_ps = fastconvAH(logicalspike , [0; PS]', size(logicalspike,1), size(logicalspike,2) );    
    lcif = lcif_mu + lcif_kx + lcif_ps;
    
    glm_ratepersec  = exp(lcif);
    glm_rateperbin  = Basepars.dt* glm_ratepersec;
    
    
    
    spikerate_bin    = size(find(logicalspike(:))) /  size(logicalspike(:));      
	model_null0      = spikerate_bin * ones(1, params.bins);
	model_uop0        = (1/params.trials) * sum(logicalspike,1);
    model_null       = repmat(model_null0, params.trials, 1);
    model_uop         = repmat(model_uop0, params.trials, 1);
    null_logprob     = sum(eval_rasterlogprob(logicalspike, model_null, 'binary', 'conditioned'));
    uop_logprob       = sum(eval_rasterlogprob(logicalspike, model_uop, 'binary', 'conditioned'));
    % Check computations are correct % 
    %null_logprob      = sum(eval_rasterlogprob(logicalspike, model_null0, 'notbinary', 'unconditioned'));
    %uop_logprob       = sum(eval_rasterlogprob(logicalspike, model_uop0, 'notbinary', 'unconditioned'));    
	uop_bits          = uop_logprob - null_logprob;
    uop_bits_perspike = uop_bits / (sum(model_null0));
    uop_bits_persecond   = uop_bits / params.testdur_seconds;
    
   
    
    
    [raster_logprob_bin] = eval_rasterlogprob( logicalspike, glm_rateperbin,  'binary', 'conditioned') ;
    glm_logprob       = sum(raster_logprob_bin);
	glm_bits          = glm_logprob - null_logprob;
    glm_bits_perspike = glm_bits / (sum(model_null0));
    glm_bits_perbin   = glm_bits / params.bins;
    glm_bits_persecond   = glm_bits / params.testdur_seconds;
    
    
    model_comparison.logprob_glm_raw(i_test, index)      = glm_logprob;
    model_comparison.logprob_uop_raw(i_test, index)      = uop_logprob;
    
    model_comparison.logprob_glm_bpspike(i_test, index)  = glm_bits_perspike;
    model_comparison.logprob_uop_bpspike(i_test, index)  = uop_bits_perspike;
    
    model_comparison.logprob_glm_bpsec(i_test, index)    = glm_bits_persecond;
    model_comparison.logprob_uop_bpsec(i_test, index)    = uop_bits_persecond;

end



    end
    model_comparison.computationtime = datestr(clock);
    model_comparison.computationmfile = mfilename('fullpath');
    
    if shortlist
        eval(sprintf('save %s/%s.mat model_comparison', savecomparisonsdir,exp_nm));
    else
        eval(sprintf('save %s/longlist_%s.mat model_comparison', savecomparisonsdir,exp_nm));
    end
    full{i_exp} = model_comparison;
    
end
fullmodelcomparison = full;

if shortlist
        eval(sprintf('save %s/%s_alldsets.mat fullmodelcomparison', savecomparisonsdir, mod_version));
else
        eval(sprintf('save %s/%s_longlist_alldsets.mat fullmodelcomparison', savecomparisonsdir, mod_version));
end


%% Some Awesome Plotting  dependent only on structure fullmodelcomparison
%
zzz = fullmodelcomparison; 
clf
subplot(3,1,1);
axis off
set(gca, 'fontsize', 12)
c = 0;
c = c+1;
text(-.1, 1-0.1*c,sprintf('Purpose is to: %s', zzz{1}.purpose ));
c = c+1;
text(-.1, 1-0.1*c,sprintf('UOP: Unconditioned Optimal Prediction, ie. the mean rate of the raster'));
c = c+1;
text(-.1, 1-0.1*c,sprintf('"Bits" refers to the (LogProb of Model - LogProb of Uniform Firing Rate)'));
c = c+1;
text(-.1, 1-0.1*c,sprintf('Each color is a different retinal prep'));
c = c+1;
text(-.1, 1-0.1*c,sprintf('Symbol * is for OFF Parasol, Symbol . is for ON Parasol'));
c = c+1;
text(-.1, 1-0.1*c,sprintf('If bits are negative, consider model fail , plotted as zero'));
c = c+1;
text(-.1, 1-0.1*c,sprintf('Computation Date %s',zzz{1}.computationtime));
c = c+1;
text(-.1, 1-0.1*c,sprintf('Plot Date %s',datestr(clock)));
c = c+2;
text(-.1, 1-0.1*c,sprintf('Computation mfile: %s',zzz{1}.computationmfile ) );

MS = 16;
for i_measure = 1:2
        amax= [];
        amin = [];
        for i_exp = 1:4
            aaa = zzz{i_exp};
            if i_measure == 1,subplot(3,2,[3 5]); measure = (aaa.logprob_glm_bpspike);    end
            if i_measure == 2,subplot(3,2,[4 6]); measure = (aaa.logprob_glm_bpspike ./ aaa.logprob_uop_bpspike);end
            

            
            ONP =  find(aaa.ONPar);
            OFFP = find(aaa.OFFPar);
            
            modelfails = find(measure<=0);
            measure(modelfails) = 0 ;
            amax = [amax max(measure(:,:))];            
            newmin = min(min(min(measure(:,ONP))) , min(min(measure(:,OFFP))));
            amin = [amin newmin];
            
            switch i_exp
                case 1
                    plot(  measure(1,OFFP),  measure(2,OFFP)  ,'r*', 'markersize' , MS ); 
                    if i_measure == 1
                        subplot(3,2,[3 5]); set(gca, 'fontsize', 12); 
                        xlabel(aaa.comparenames{1}); ylabel(aaa.comparenames{2}); title('GLM Perfomance: Bits Per Spike')  ; 
                    end
                    if i_measure == 2
                        subplot(3,2,[4 6]); set(gca, 'fontsize', 12); 
                        xlabel(aaa.comparenames{1}); ylabel(aaa.comparenames{2}); title('GLMperformance / UOPperformance'); 
                    end
                    hold on
                    plot(  measure(1,ONP),   measure(2,ONP)  ,'r.', 'markersize' , MS ); 
                case 2
                    plot(  measure(1,OFFP),  measure(2,OFFP)  ,'g*', 'markersize' , MS ); 
                    plot(  measure(1,ONP),  measure(2,ONP)  ,'g.', 'markersize' , MS ); 
                case 3
                    plot(  measure(1,OFFP),  measure(2,OFFP)  ,'b*', 'markersize' , MS ); 
                    plot(  measure(1,ONP),  measure(2,ONP)  ,'b.', 'markersize' , MS ); 
                case 4
                    plot(  measure(1,OFFP),  measure(2,OFFP)  ,'k*', 'markersize' , MS ); 
                    plot(  measure(1,ONP),  measure(2,ONP)  ,'k.', 'markersize' , MS ); 
            end
            
        end
        zmax = max(amax);
        zmin = min(amin);
        plot(linspace(zmin,zmax,100), linspace(zmin,zmax,100),'k');
        xlim([zmin zmax]) ; ylim([zmin zmax]);
end
orient landscape
eval(sprintf('print -dpdf %s/%s.pdf', zzz{1}.savedir, zzz{1}.mod_version));
%}            

    

