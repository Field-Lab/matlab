% AKHeitman 2013-03-29
% New version should be a bit more modular

% AKHeitman 2014-03-23
% Effectively Done as of 2014-03-23
% Making Robust Model Comparison 
% Decided to fix on  (LogProb model - LogProb null) / spikes
% Normalized by (LogProb UOP - LogProb null) / spikes



%%%% Set up directories %%%%
%clear; close all; clear all
%GLMdir ='/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/GLM' ;
%cellselectiontype = 'shortlist';
%exptests = [1 2 3 4];
%fixed.fitversion = 'Fit';fixed.map_type = 'mapPRJ';
%modulation  = 'model_parameters';
%mod_version = 'WN_ONFFvsLin_8pix_ID_8pix'; purpose = 'Examine effect of ONOff Hard rectification channels in WN, no cone model';
%comp_params{1}.name = 'single linear channel';
%comp_params{2}.name = 'hard rect ON and Off channel';
%comp_params{1}.filtertype = 'fixedSP';
%comp_params{2}.filtertype = 'OnOff_hardrect_fixedSP_STA';
%comp_params{1}.cmodel    ='8pix_Identity_8pix';
%comp_params{2}.cmodel    ='8pix_Identity_8pix';
%fixed.fit_type = 'BW'; fixed.test_type = 'BW';






%%% MODULATION AND MOD_VERSION %%%
% define  changed and what stays fixed
% modulation is either, "fit_and_test" or "model_parameters" or "other"
% "fit_and_test" referes to NSEM or WN for fitting and testing
% "model_parameters"  should be flexible for any model change




%%%% FIT AND TEST CHANGES %%%%%
%{





%modulation = 'fit_and_test';
mod_version = 'WNWN_vs_NSEMNSEM'; purpose = 'Show WN GLM fits better than NSEM GLM fits'; 
comp_params{1}.name = 'fit/test by WN'; 
comp_params{2}.name = 'fit/test by NSEM';

%mod_version = 'WNWN_vs_WNNSEM'; purpose = 'Show WN fits do not generalize to NSEM'; 
%comp_params{1}.name = 'fit/test by WN'; 
%comp_params{2}.name = 'fit with WN, test with NSEM';

%mod_version = 'NSEMWN_vs_NSEMNSEM'; purpose = 'Show NSEM fits do not generalize to NSEM'; 
%comp_params{1}.name = 'fit with NSEM, test with WN';
%comp_params{2}.name = 'fit/test by NSEM'; 
modulation  = 'model_parameters';
mod_version = 'WN_conemodel_8pix_1e4_8pix'; purpose = 'Show Cone model has minimal effect for WN fits';
comp_params{1}.name = 'fit without cone model';
comp_params{2}.name = 'fit with cone model 8p1e48p';
comp_params{1}.filtertype = 'fixedSP';
comp_params{2}.filtertype = 'fixedSP';
comp_params{1}.cmodel    ='8pix_Identity_8pix';
comp_params{2}.cmodel    ='8pix_Model1_1e4_8pix';
fixed.fit_type = 'BW'; fixed.test_type = 'BW';
%}

%%% Write up savedir creat saving structure %%%

function fullmodelcomparison = comparefits_LOGPROB_robust(mod, comp_params, fixed, exptests , cellselectiontype, GLMdir,plots)

purpose = mod.purpose;
mod_version = mod.mod_version;
modulation = mod.modulation

modelcompare_basedir = sprintf('%s/Model_Compare',GLMdir)
savecomparisonsdir = sprintf('%s/%s/%s', modelcompare_basedir , modulation, mod_version);
if ~exist(savecomparisonsdir, 'dir'), mkdir(savecomparisonsdir); end
full = cell( length(exptests) ,1);
i_exp = 1; i_test = 1;
%%
for i_exp = exptests 
    %% Setup model_comparison structure  one per dataset
    
    [exp_nm,cells,expname]  = cell_list( i_exp, cellselectiontype);
    head_cellID = [cells{:}];
    
    clear model_comparison
    model_comparison.exp_nm            = exp_nm;
    model_comparison.exp_nm            = exp_nm;
    model_comparison.cell_ids          = head_cellID;   
    model_comparison.comparison_type   = modulation;
    model_comparison.mod_version       = mod_version;
    model_comparison.purpose           = purpose;
    model_comparison.comparenames      = { comp_params{1}.name , comp_params{2}.name };
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
        map_type = fixed.map_type ;   
        fittype = sprintf('%s_%s/%s_%s' , fit_type, map_type, fit_type, fixed.fitversion);  

        [StimulusPars DirPars datarun_slv datarun_mas] = Directories_Params_v20(exp_nm, fit_type, 'mapPRJ');
        if strcmp(test_type , 'NSEM')
            SPars = StimulusPars.NSEM;
        elseif strcmp(test_type, 'BW')
            SPars = StimulusPars.BW;
        end
        orgspikesdir = sprintf('/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/GLM/%s/%s_mapPRJ/BlockSpikesandRasters', exp_nm, test_type);  
        stimulidir = '/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Stimuli';
        [blockedmoviecell,  novelmoviestats, origmatfile] = loadmoviematfile(exp_nm , test_type, cmodel, 'testmovie')  ; 
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
        cid = head_cellID(1);

%% Individual Cell Loops . try to make this a functional unit for evaluation w/n GLM fitting
for cid = head_cellID
    
    %%%% LOAD BASEPARS %%%%%%%%%%%%%%%
    %%%% LOAD ORGANIZED SPIKES %%%%%%%
    index = find( head_cellID == cid);
	if ~isempty(find(datarun_mas.cell_types{1}.cell_ids == cid))
            celltype = 'ONPar';
            model_comparison.ONPar(index) = 1;
            
    end
	if ~isempty(find(datarun_mas.cell_types{2}.cell_ids == cid))
            celltype = 'OFFPar';
            model_comparison.OFFPar(index) = 1;
    end
    params.sname = sprintf('%s_%d', celltype, cid);
    display(sprintf('Will evaluate xval performance for %s', params.sname))
    load(sprintf('%s/%s/%s/%s/%s.mat' , GLMdir, exp_nm,fittype, fitparams , params.sname));
    eval(sprintf('load %s/organizedspikes_%s.mat', orgspikesdir , params.sname));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% EVALUATE PERFORMANCE %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xvalperformance = eval_xvalperformance(Basepars, SPars, organizedspikes,testmovie, novelmoviestats.mu_avgIperpix);    
    
    model_comparison.logprob_null_raw(i_test,index)      =  xvalperformance.logprob_null_raw ;
    model_comparison.logprob_uop_raw(i_test,index)       =  xvalperformance.logprob_uop_raw ;
    model_comparison.logprob_glm_raw(i_test,index)       =  xvalperformance.logprob_glm_raw;
    model_comparison.logprob_uop_bpspike(i_test,index)   =  xvalperformance.logprob_uop_bpspike; 
    model_comparison.logprob_glm_bpspike(i_test,index)   =  xvalperformance.logprob_glm_bpspike;
    model_comparison.logprob_uop_bpsec(i_test,index)     =  xvalperformance.logprob_uop_bpsec;
    model_comparison.logprob_glm_bpsec(i_test,index)     =  xvalperformance.logprob_glm_bpsec ;
    model_comparison.glm_normedbits(i_test,index)        =  xvalperformance.glm_normedbits;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end



    end
    model_comparison.computationtime = datestr(clock);
    model_comparison.computationmfile = mfilename('fullpath');
    
    
    eval(sprintf('save %s/%s_%s.mat model_comparison', savecomparisonsdir,cellselectiontype , exp_nm));

    full{i_exp} = model_comparison;
    
end
fullmodelcomparison = full;

eval(sprintf('save %s/%s_crossprep.mat fullmodelcomparison', savecomparisonsdir, cellselectiontype));


%% Some Awesome Plotting  dependent only on structure fullmodelcomparison

if plots
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
        for i_exp = exptests 
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
eval(sprintf('print -dpdf %s/%s_%s.pdf', zzz{1}.savedir,cellselectiontype, zzz{1}.mod_version));
end




end
            

    

