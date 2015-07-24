% akheitman 2014-03-08 .. finished 2014-03-09
% Rater_Analysis_logicalspike
% clean up the COSYNE scramble
% Good Solid Quantification of Rasters  .. no model yet

% will ty to bin spikes not just logical.. can have 2 3 4 etc. spikes

clear; close all; clear all
GLMdir ='/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/GLM' ;
figuregeneratingfilename = mfilename('fullpath');
shortlist = false;
savedir = sprintf('/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/Rasters_WNandNSEM/logprob_predictability');
if ~exist(savedir), mkdir(savedir); end

for bpf = [.25 .5 1 2 5 10];
tstim = .0083275;
bindur = tstim / bpf; 

CrossPrep_RasterAnal = cell(1,4);    
%%
for i_exp = 1:4
    %% 
    if i_exp == 1, exp_nm = '2012-08-09-3' ; expname = 'expA'; head_cellID =   [1471  841 3676 5086 5161 1426 1772 2101 1276]; end
	if i_exp == 2, exp_nm = '2012-09-27-3';  expname = 'expB'; head_cellID =     [1 31 301 1201 1726 91 1909 2360 6858];  end
	if i_exp == 3, exp_nm = '2013-08-19-6';  expname = 'expC'; head_cellID = [737 1328 1341 2959  2824 3167 3996 5660 6799 5447]; end
	if i_exp == 4, exp_nm = '2013-10-10-0';  expname = 'expD'; head_cellID = [346 1233 3137 7036 5042 5418 32 768 2778 4354 5866];end
    [StimulusPars DirPars datarun_NSEM datarun_mas] = Directories_Params_v20(exp_nm, 'NSEM', 'mapPRJ');
    [StimulusPars DirPars datarun_BW   datarun_mas] = Directories_Params_v20(exp_nm, 'BW', 'mapPRJ');

    
    if ~shortlist 
        clear head_cellID; 
        head_cellIDBW = [datarun_BW.cell_types{2}.cell_ids datarun_BW.cell_types{1}.cell_ids];
        head_cellIDNSEM = [datarun_NSEM.cell_types{2}.cell_ids datarun_NSEM.cell_types{1}.cell_ids];
        
        head_cellID = intersect(head_cellIDBW , head_cellIDNSEM)
    end

  
    
    rasteranal.exp_nm    = exp_nm;
    rasteranal.binsperframe   = bpf;
    rasteranal.binduration_secs = bindur;
    rasteranal.cell_ids  = head_cellID;
    rasteranal.note   = 'reflect the timing of model evaluation.. might not be full raster duration (in neither trials nor seconds) '
    rasteranal.note2  = 'first row is WN numbers ,  second row is NSEM numbers';
    rasteranal.mr_logprobs                 = zeros(2,length(head_cellID) ); 
    rasteranal.mr_logprob_perbin           = zeros(2,length(head_cellID) ); 
    rasteranal.null_logprob                = zeros(2,length(head_cellID) ); 
    rasteranal.null_logprob_perbin         = zeros(2,length(head_cellID) ); 
    rasteranal.mr_bits                     = zeros(2,length(head_cellID) ); 
    rasteranal.mr_bits_perbin              = zeros(2,length(head_cellID) ); 
    rasteranal.mr_bits_perspike            = zeros(2,length(head_cellID) ); 
    rasteranal.firingrate_persec           = zeros(2,length(head_cellID) ); 
    rasteranal.firingrate_perbin           = zeros(2,length(head_cellID) ); 
    rasteranal.ONPar= zeros(1,length(head_cellID) ); 
    rasteranal.OFFPar = zeros(1,length(head_cellID) ); 
    cid = head_cellID(1);
    i_stimtype = 1;
    
%%
    for cid = head_cellID 
              
        for i_stimtype = 1:2  
            
        %%
        if i_stimtype == 1
            stim_type = 'BW';
            SPars = StimulusPars.BW;
        elseif i_stimtype == 2
            stim_type = 'NSEM';
            SPars = StimulusPars.NSEM;           
        end
        evalblocks = SPars.evalmodel_Blocks;
        orgspikesdir = sprintf('/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects/GLM/%s/%s_mapPRJ/BlockSpikesandRasters', exp_nm, stim_type); 
        
        index = find( head_cellID == cid);
        if ~isempty(find(datarun_mas.cell_types{1}.cell_ids == cid))
            celltype = 'ONPar';
            rasteranal.ONPar(index) = 1;
        end
        if ~isempty(find(datarun_mas.cell_types{2}.cell_ids == cid))
            celltype = 'OFFPar';
            rasteranal.OFFPar(index) = 1;
        end    
        sname = sprintf('%s_%d', celltype, cid);
        eval(sprintf('load %s/organizedspikes_%s.mat', orgspikesdir , sname));
        
        %Assumes the bin size small enough to only have 1 spike per bin
        testbins = bpf *length(SPars.testframes);   
        trials = length(evalblocks);       
        logicalspike = zeros( length(evalblocks) , testbins) ; 
        
        for i_blk = 1 : length(evalblocks)
            blknum = evalblocks(i_blk);
            sptimes = organizedspikes.block.t_sp_withinblock{blknum} - SPars.testseconds(1);
            sptimes = sptimes(find(sptimes > 0 ) );
            for i_sp = 1:length(sptimes)
                spt = sptimes(i_sp);
                binnumber = ceil(spt / bindur );
                logicalspike( i_blk, binnumber )  =  logicalspike( i_blk,binnumber ) + 1;
            end
        end
        
        if i_stimtype == 1
            BWrast = logicalspike;
        elseif i_stimtype == 2
            NSEMrast = logicalspike;
        end
        %%
        spikerate_bin    = size(find(logicalspike(:))) /  size(logicalspike(:));      
        model_null       = spikerate_bin * ones(1, testbins);
        model_mr         = (1/trials) * sum(logicalspike,1);       
        null_logprob     = sum(eval_rasterlogprob(logicalspike, model_null, 'notbinary', 'unconditioned'));
        mr_logprob       = sum(eval_rasterlogprob(logicalspike, model_mr, 'notbinary', 'unconditioned'));
        mr_bits          = mr_logprob - null_logprob;
        mr_bits_perspike = mr_bits / (sum(model_null));
        mr_bits_perbin   = mr_bits / testbins;



        rasteranal.mr_logprobs(i_stimtype, index)                  = mr_logprob;
        rasteranal.mr_logprob_perbin(i_stimtype, index)            = mr_logprob / testbins;
        rasteranal.mr_logprob_perspike(i_stimtype, index)          = mr_logprob / (sum(model_null)) ;
        rasteranal.null_logprob(i_stimtype, index)                 = null_logprob;
        rasteranal.null_logprob_perbin(i_stimtype, index)          = null_logprob / testbins;
        rasteranal.null_logprob_perspike(i_stimtype, index)        = null_logprob / (sum(model_null));
        rasteranal.mr_bits(i_stimtype, index)                      = mr_bits; 
        rasteranal.mr_bits_perbin(i_stimtype, index)               = mr_bits_perbin;
        
        rasteranal.mr_bits_perspike(i_stimtype, index)             = mr_bits_perspike;
        
        
        
        rasteranal.firingrate_persec(i_stimtype, index)            = spikerate_bin / bindur;
        rasteranal.firingrate_perbin(i_stimtype, index)            = spikerate_bin;
      
        
        end
    end
    CrossPrep_RasterAnal{i_exp} = rasteranal;
end
if shortlist
    eval(sprintf('save %s/slist_RasterAnal_binsperframe_%d.mat CrossPrep_RasterAnal' , savedir, bpf));
end

%%
clf;
subplot(3,1,1);
axis off
set(gca, 'fontsize', 12)
c = 0;
c = c+1;
text(-.1, 1-0.1*c,sprintf('Bin duration of %1.4e millisecs', 1000*bindur));
c = c+1;
text(-.1, 1-0.1*c,sprintf('UOP: Unconditioned Optimal Prediction, ie. the mean rate of the raster, binned at the bin duratiion'));
c = c+1;
text(-.1, 1-0.1*c,sprintf('Quantitative Analaysis of the assertion that Natural Scenes rasters are  "tighter" and less noisy than WN rasters'));
c = c+1;
text(-.1, 1-0.1*c,sprintf('"Bits" refers to the (LogProb of Model - LogProb of Steady Firing Rate)'));
c = c+1;
text(-.1, 1-0.1*c,sprintf('Each Color is a different retinal prep'));
c = c+1;
text(-.1, 1-0.1*c,sprintf('%s', datestr(clock)));
c = c+1;
text(-.1, 1-0.1*c,sprintf('%s',figuregeneratingfilename ) );

if shortlist
   MS = 16;
else
   MS = 8;
end
for i_measure = 1:4
        if i_measure == 1,subplot(3,2,3); set(gca, 'fontsize', 12); xlabel('WhiteNoise'); ylabel('Natural Scenes'); title('UOP Perfomance: LogProb PerBin')  ; end
        if i_measure == 2,subplot(3,2,4); set(gca, 'fontsize', 12); xlabel('WhiteNoise'); ylabel('Natural Scenes'); title('UOP Perfomance: LogProb PerSpike'); end
        if i_measure == 3,subplot(3,2,5); set(gca, 'fontsize', 12); xlabel('WhiteNoise'); ylabel('Natural Scenes'); title('UOP Perfomance: Bits PerBin'); end
        if i_measure == 4,subplot(3,2,6); set(gca, 'fontsize', 12); xlabel('WhiteNoise'); ylabel('Natural Scenes'); title('UOP Perfomance: Bits PerSpike');end 
        amax= [];
        amin = [];
        for i_exp = 1:4
            rasteranal = CrossPrep_RasterAnal{i_exp};
            if i_measure == 1,subplot(3,2,3); measure = rasteranal.mr_logprob_perbin; end
            if i_measure == 2,subplot(3,2,4); measure = rasteranal.mr_logprob_perspike;end
            if i_measure == 3,subplot(3,2,5); measure = rasteranal.mr_bits_perbin;end
            if i_measure == 4,subplot(3,2,6); measure = rasteranal.mr_bits_perspike;end
            amax = [amax max(measure(:))];
            amin = [amin min(measure(:))];
            
            ONP =  find(rasteranal.ONPar);
            OFFP = find(rasteranal.OFFPar);
            switch i_exp
            case 1
                hold on;
                plot(  measure(1,OFFP),  measure(2,OFFP)  ,'r.', 'markersize' , MS ); 
                plot(  measure(1,ONP),   measure(2,ONP)  ,'r.', 'markersize' , MS ); 
            case 2
                plot(  measure(1,OFFP),  measure(2,OFFP)  ,'g.', 'markersize' , MS ); 
                plot(  measure(1,ONP),  measure(2,ONP)  ,'g.', 'markersize' , MS ); 
            case 3
                plot(  measure(1,OFFP),  measure(2,OFFP)  ,'b.', 'markersize' , MS ); 
                plot(  measure(1,ONP),  measure(2,ONP)  ,'b.', 'markersize' , MS ); 
            case 4
                plot(  measure(1,OFFP),  measure(2,OFFP)  ,'k.', 'markersize' , MS ); 
                plot(  measure(1,ONP),  measure(2,ONP)  ,'k.', 'markersize' , MS ); 
            end
        end
        zmax = max(amax);
        zmin = min(amin);
        plot(linspace(zmin,zmax,100), linspace(zmin,zmax,100),'k');
        xlim([zmin zmax]) ; ylim([zmin zmax]);
end
orient tall
if shortlist
    filename = sprintf('Raster_WNvsNSEM_OptimalPrediction_Bindur_%1.2e_msecs_subset', 1000*bindur);
else
    filename = sprintf('Raster_WNvsNSEM_OptimalPrediction_Bindur_%1.2e_msecs', 1000*bindur);
end
eval(sprintf('print -dpdf %s/%s.pdf', savedir, filename))
end
