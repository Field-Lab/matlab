% AKHEITMAN 2014-11-11
% New Version Includes More TIme Scales 
% Wrapper for Raster computations
% Used to decide which cells to keep

clear; close all; clc
BD = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection));  % allcells structire
expanded_metrics = true;
rastmetdir = sprintf('%s/raster_performance', BD.Cell_Selection);
if ~exist(rastmetdir, 'dir'), mkdir(rastmetdir); end

raster_scores = allcells;

DirPars.rast_dir = BD.BlockedSpikes;
binscales    = [1 2 5 10 25 50 100];

params.bindur     = .00083275;
params.bins_per_frame = 10;
metricparams.smoothbins    = binscales;
metricparams.bindur = params.bindur;
if exist('expanded_metrics', 'var') && expanded_metrics
    metricparams.pair_numbers      = 40;
end
    
i_exp = 1; i_stimtype = 1; i_celltype = 1; i_cell = 1;

individual_plots = true;
%exps = 1:2
exps = 3:4
%%

for i_exp = 1:4
    raster_scores{i_exp}.scores_stimtype_ind = 'One is White Noise, Two is NSEM';
    raster_scores{i_exp}.scores_celltype_ind = 'One is On Parasols, Two is Off Parasols';
    raster_scores{i_exp}.expand_binscales = binscales;
    raster_scores{i_exp}.expand_timescales = params.bindur*binscales;
end


%% LOOP THROUGH AND COMPUTE RASTER SCORES
for i_exp = exps
    for i_stimtype = 1:2    
        for i_celltype = 1:2
            %%
            exp_nm  = allcells{i_exp}.exp_nm;
            expname = allcells{i_exp}.expname;
            if i_stimtype == 1, stimtype = 'WN';   end
            if i_stimtype == 2, stimtype = 'NSEM'; end
            inputs.exp_nm       = exp_nm; 
            inputs.map_type     = 'mapPRJ';
            inputs.stim_type    = stimtype;
            [StimulusPars DirPars] = Directories_Params_v23(exp_nm, inputs.stim_type, inputs.map_type);
             DirPars.organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', inputs);
             SPars = StimulusPars.slv; 
             params.evalblocks = SPars.TestBlocks;
             params.bins = params.bins_per_frame *length(SPars.testframes); 
            if i_celltype == 1; cells = allcells{i_exp}.ONP;  celltype = 'ONPar'; end
            if i_celltype == 2; cells = allcells{i_exp}.OFFP; celltype = 'OFFPar'; end
            
            
            
            
            % Set up structure
            if exist('expanded_metrics','var') && expanded_metrics
                for i_t = 1:length(binscales)
                    raster_scores{i_exp}.expand{i_t}.stim_type{i_stimtype}.celltype{i_celltype}.temporal_variance_normed = zeros(1,length(cells));
                    raster_scores{i_exp}.expand{i_t}.stim_type{i_stimtype}.celltype{i_celltype}.fractional_error = zeros(1,length(cells));
                    raster_scores{i_exp}.expand{i_t}.stim_type{i_stimtype}.celltype{i_celltype}.viktordist_perspike = zeros(1,length(cells));
                    raster_scores{i_exp}.expand{i_t}.stim_type{i_stimtype}.celltype{i_celltype}.bps_mean = zeros(1,length(cells));
                    raster_scores{i_exp}.expand{i_t}.stim_type{i_stimtype}.celltype{i_celltype}.bps_max = zeros(1,length(cells));
                    raster_scores{i_exp}.expand{i_t}.stim_type{i_stimtype}.celltype{i_celltype}.bps_std = zeros(1,length(cells));
                    raster_scores{i_exp}.expand{i_t}.stim_type{i_stimtype}.celltype{i_celltype}.logprob_persec_mean = zeros(1,length(cells));
                    raster_scores{i_exp}.expand{i_t}.stim_type{i_stimtype}.celltype{i_celltype}.logprob_persec_max = zeros(1,length(cells));
                    raster_scores{i_exp}.expand{i_t}.stim_type{i_stimtype}.celltype{i_celltype}.logprob_persec_std = zeros(1,length(cells));
                end
            else
                raster_scores{i_exp}.scores.stim_type{i_stimtype}.celltype{i_celltype}.spikespersec      = zeros(1,length(cells)); 
                raster_scores{i_exp}.scores.stim_type{i_stimtype}.celltype{i_celltype}.trial_variance    = zeros(1,length(cells));
                raster_scores{i_exp}.scores.stim_type{i_stimtype}.celltype{i_celltype}.temporal_variance = zeros(1,length(cells));
                raster_scores{i_exp}.scores.stim_type{i_stimtype}.celltype{i_celltype}.bits_per_spike    = zeros(1,length(cells));
            end
            
            
            
            
            %%
            for i_cell = 1:length(cells)
                %% 
                % LOAD CELLS INFORMATION AND GENERATE RASTER %
                cid = cells(i_cell);
                cell_savename = sprintf('%s_%d', celltype,cid);
                display(sprintf('Working on %s  Cell: %s', exp_nm,cell_savename));
                eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', DirPars.organizedspikesdir, cell_savename));
                %%%  Fill in the logical spikes with ~msec bins %%%
                logicalspike = zeros( length(params.evalblocks) , params.bins) ;         
                for i_blk = 1 : length(params.evalblocks)
                    blknum = params.evalblocks(i_blk);
                    sptimes = organizedspikes.block.t_sp_withinblock{blknum} - SPars.fittest_skipseconds;
                    sptimes = sptimes(find(sptimes > 0 ) );
                    for i_sp = 1:length(sptimes)
                        spt = sptimes(i_sp);
                        binnumber = ceil(spt / params.bindur );
                        logicalspike( i_blk, binnumber )  =  logicalspike( i_blk,binnumber ) + 1;
                    end
                end
                
                % DO RASTER COMPUTATIONS %
                if exist('expanded_metrics','var') && expanded_metrics
                    ind_raster_scores = raster_metrics_expanded(logicalspike,metricparams);
                else
                    tic
                    ind_raster_scores = raster_metrics(logicalspike,metricparams);
                    toc
                end
                irs = ind_raster_scores;
                
                
                
                if exist('expanded_metrics','var') && expanded_metrics
                    for i_t = 1:length(binscales)
                        
                        raster_scores{i_exp}.expand{i_t}.stim_type{i_stimtype}.celltype{i_celltype}.sigma_bin = irs{i_t}.param_sigma_bin;
                        raster_scores{i_exp}.expand{i_t}.stim_type{i_stimtype}.celltype{i_celltype}.sigma_sec = irs{i_t}.param_sigma_bin *params.bindur ;
                        raster_scores{i_exp}.expand{i_t}.stim_type{i_stimtype}.celltype{i_celltype}.temporal_variance_normed(i_cell) = irs{i_t}.temporal_variance_normed;
                        raster_scores{i_exp}.expand{i_t}.stim_type{i_stimtype}.celltype{i_celltype}.fractional_error(i_cell) = irs{i_t}.fractional_error;
                        raster_scores{i_exp}.expand{i_t}.stim_type{i_stimtype}.celltype{i_celltype}.viktordist_perspike(i_cell) = irs{i_t}.viktordist_perspike;
                        raster_scores{i_exp}.expand{i_t}.stim_type{i_stimtype}.celltype{i_celltype}.bps_mean(i_cell) = irs{i_t}.bps_mean; 
                        raster_scores{i_exp}.expand{i_t}.stim_type{i_stimtype}.celltype{i_celltype}.bps_max(i_cell) = irs{i_t}.bps_max;
                        raster_scores{i_exp}.expand{i_t}.stim_type{i_stimtype}.celltype{i_celltype}.bps_std(i_cell) = irs{i_t}.bps_std;
                        raster_scores{i_exp}.expand{i_t}.stim_type{i_stimtype}.celltype{i_celltype}.logprob_persec_mean(i_cell) = irs{i_t}.logprobpersec_mean; 
                        raster_scores{i_exp}.expand{i_t}.stim_type{i_stimtype}.celltype{i_celltype}.logprob_persec_max(i_cell) = irs{i_t}.logprobpersec_max;
                        raster_scores{i_exp}.expand{i_t}.stim_type{i_stimtype}.celltype{i_celltype}.logprob_persec_std(i_cell) = irs{i_t}.logprobpersec_std;

                        
                    end
                else
                    raster_scores{i_exp}.scores.stim_type{i_stimtype}.celltype{i_celltype}.spikespersec(i_cell)      = irs.spikespersec;
                    raster_scores{i_exp}.scores.stim_type{i_stimtype}.celltype{i_celltype}.trial_variance(i_cell)    = irs.trial_variance;
                    raster_scores{i_exp}.scores.stim_type{i_stimtype}.celltype{i_celltype}.temporal_variance(i_cell) = irs.temporal_variance;
                    raster_scores{i_exp}.scores.stim_type{i_stimtype}.celltype{i_celltype}.bits_per_spike(i_cell)    = irs.bits_per_spike;
                end

                
            end     
        end
    end
    
    
    if i_exp == 1,expanded_raster_scoresA = raster_scores{i_exp}; eval(sprintf('save %s/expanded_raster_scoresA.mat expanded_raster_scoresA',rastmetdir)); end
    if i_exp == 2,expanded_raster_scoresB = raster_scores{i_exp}; eval(sprintf('save %s/expanded_raster_scoresB.mat expanded_raster_scoresB',rastmetdir)); end
    if i_exp == 3,expanded_raster_scoresC = raster_scores{i_exp}; eval(sprintf('save %s/expanded_raster_scoresC.mat expanded_raster_scoresC',rastmetdir)); end
    if i_exp == 4,expanded_raster_scoresD = raster_scores{i_exp}; eval(sprintf('save %s/expanded_raster_scoresD.mat expanded_raster_scoresD',rastmetdir)); end
end

%{
if exist('expanded_metrics','var') && expanded_metrics
    expanded_raster_scores = raster_scores;
	eval(sprintf('save %s/expanded_raster_scores.mat expanded_raster_scores',rastmetdir))
else
	eval(sprintf('save %s/raster_scores.mat raster_scores',rastmetdir))
end
%}
%% SCATTER PLOTS OF THE RASTER SCORES
%{
eval(sprintf('load %s/expanded_raster_scoresA.mat', rastmetdir));
eval(sprintf('load %s/expanded_raster_scoresB.mat', rastmetdir));
eval(sprintf('load %s/expanded_raster_scoresC.mat', rastmetdir));
eval(sprintf('load %s/expanded_raster_scoresD.mat', rastmetdir));

close all
expts = 1:4
i_metric = 6 ; i_time = 1; i_exp = 1; i_type = 1;
%%
for i_metric = 1:6
    for i_time = 1:3
        clf; hold on; set(gca,'fontsize', 12)
        if i_metric == 1, title_base = 'Temporal Variance Normed';  pdf_title = 'TemporalVariance';     end
        if i_metric == 2, title_base = 'Bits Per Spike Mean';       pdf_title = 'BPS_Mean';             end
        if i_metric == 3, title_base = 'Bits Per Spike Max';        pdf_title = 'BPS_Max';              end
        if i_metric == 4, title_base = 'Bits Per Spike Std';        pdf_title = 'BPS_Std';              end
        if i_metric == 5, title_base = 'Fractional Error Capped at 30';          pdf_title = 'pair_distance_vector'; end
        if i_metric == 6, title_base = 'Viktor Distance Per Spike'; pdf_title = 'pair_distance_Viktor'; end
        
        if i_time ==1 , bins = metricparams.smoothbins; end
        if i_time ==2 , bins = metricparams.slowsmoothbins; end
        if i_time ==3 , bins = metricparams.slowestsmoothbins; end
        timeadded = sprintf('TimeScale %d msecs', round(1000*params.bindur*bins) );
        pdftime = sprintf('timescale_%dBins', bins);
        
        title(sprintf('%s %s', title_base, timeadded));
        

        xlabel('WN values');
        ylabel('NSEM values');


        MS = 12;
        max_x = 0;
        max_y = 0;
        for i_exp = expts
            if i_exp == 1; basecolor = 'r'; fullscores = expanded_raster_scoresA.expand{i_time}; end
            if i_exp == 2; basecolor = 'g'; fullscores = expanded_raster_scoresB.expand{i_time}; end
            if i_exp == 3; basecolor = 'b'; fullscores = expanded_raster_scoresC.expand{i_time}; end
            if i_exp == 4; basecolor = 'c'; fullscores = expanded_raster_scoresD.expand{i_time}; end

            for i_type = 1:2
                if i_type == 1, marktype  = '.'; end
                if i_type == 2, marktype  = '*'; end

                wn_struct   = fullscores.stim_type{1}.celltype{i_type};
                nsem_struct = fullscores.stim_type{2}.celltype{i_type};

                if i_metric == 1, wn_vec = sqrt(wn_struct.temporal_variance_normed); end
                if i_metric == 2, wn_vec = wn_struct.bps_mean;    end
                if i_metric == 3, wn_vec = wn_struct.bps_max; end
                if i_metric == 4, wn_vec = wn_struct.bps_std;    end
                if i_metric == 5, wn_vec = wn_struct.fractional_error;    end
                if i_metric == 6, wn_vec = wn_struct.viktordist_perspike;    end

                if i_metric == 1, nsem_vec = sqrt(nsem_struct.temporal_variance_normed); end
                if i_metric == 2, nsem_vec = nsem_struct.bps_mean;    end
                if i_metric == 3, nsem_vec = nsem_struct.bps_max; end
                if i_metric == 4, nsem_vec = nsem_struct.bps_std;    end
                if i_metric == 5, nsem_vec = nsem_struct.fractional_error;    end
                if i_metric == 6, nsem_vec = nsem_struct.viktordist_perspike;    end
                
                
                if i_metric == 5, wn_vec(wn_vec>30) = 30; nsem_vec(nsem_vec>30)=30; end
                


                plotstring = sprintf('%s%s',basecolor,marktype);
                plot(wn_vec, nsem_vec, plotstring,'markersize',MS);             
                max_x = max(max_x, max(wn_vec));
                max_y = max(max_y, max(nsem_vec));
            end
        end
        max_val = max(max_x,max_y);
        set(gca,'xlim',[0 max_val]); set(gca,'ylim',[0 max_val]);

        plot(linspace(0,max_val,100), linspace(0,max_val,100),'k');
        orient landscape
        
        
        
        eval(sprintf('print -dpdf %s/%s_%s.pdf',rastmetdir, pdf_title,pdftime)); 
    end
end
%}

%{
close all
for i_metric = 1:4
    figure; hold on;
    if i_metric == 1, title('Spike Rate'); end
    if i_metric == 2, title('Trial Variance'); end
    if i_metric == 3, title('Temporal Variance'); end
    if i_metric == 4, title('Bits Per Spike'); end
    
    xlabel('WN values');
    ylabel('NSEM values');
    
    
    MS = 12;
    max_x = 0;
    max_y = 0;
    for i_exp = 1:4
        if i_exp == 1; basecolor = 'r'; end
        if i_exp == 2; basecolor = 'g'; end
        if i_exp == 3; basecolor = 'b'; end
        if i_exp == 4; basecolor = 'c'; end
        
        

        for i_type = 1:2
            if i_type == 1, marktype  = '.'; end
            if i_type == 2, marktype = '*' ; end
            
            wn_struct   = raster_scores{i_exp}.scores.stim_type{1}.celltype{i_type};
            nsem_struct = raster_scores{i_exp}.scores.stim_type{2}.celltype{i_type};
            
            if i_metric == 1, wn_vec = wn_struct.spikespersec;      end
            if i_metric == 2, wn_vec = wn_struct.trial_variance;    end
            if i_metric == 3, wn_vec = wn_struct.temporal_variance; end
            if i_metric == 4, wn_vec = wn_struct.bits_per_spike;    end
            
            if i_metric == 1, nsem_vec = nsem_struct.spikespersec;      end
            if i_metric == 2, nsem_vec = nsem_struct.trial_variance;    end
            if i_metric == 3, nsem_vec = nsem_struct.temporal_variance; end
            if i_metric == 4, nsem_vec = nsem_struct.bits_per_spike;    end
            

            plotstring = sprintf('%s%s',basecolor,marktype);
            plot(wn_vec, nsem_vec, plotstring,'markersize',MS);             
            max_x = max(max_x, max(wn_vec));
            max_y = max(max_y, max(nsem_vec));
        end
    end
    max_val = max(max_x,max_y);
    set(gca,'xlim',[0 max_val]); set(gca,'ylim',[0 max_val]);
    
    plot(linspace(0,max_val,100), linspace(0,max_val,100),'k');
    orient landscape
    if i_metric == 1, eval(sprintf('print -dpdf %s/rastermetric_spikerate_scatteronly.pdf',rastmetdir)); end
    if i_metric == 2, eval(sprintf('print -dpdf %s/rastermetric_trialvariance_scatteronly.pdf',rastmetdir)); end
    if i_metric == 3, eval(sprintf('print -dpdf %s/rastermetric_temporalvariance_scatteronly.pdf',rastmetdir)); end
    if i_metric == 4, eval(sprintf('print -dpdf %s/rastermetric_uop_bps_scatteronly.pdf',rastmetdir)); end
    
    
end


%% Plot WN and NSEM rasters with corresponding scores
if exist('individual_plots', 'var') && individual_plots

plotsecs = 8; MS = 6;
eval(sprintf('load %s/raster_scores.mat raster_scores',rastmetdir))
i_exp = 1; i_stimtype = 1; i_celltype = 1; i_cell = 1;
for i_exp = 1:4
    %%
    exp_nm  = allcells{i_exp}.exp_nm;
    expname = allcells{i_exp}.expname;
    
    plotdir = sprintf('%s/plot_raster_scores/%s',rastmetdir,exp_nm);
    if ~exist(plotdir,'dir'), mkdir(plotdir); end
    inputs.exp_nm       = exp_nm; 
    inputs.map_type     = 'mapPRJ';
    % Directory Loading
    for i_stimtype = 1:2
        if i_stimtype == 1, stimtype = 'WN';   end
        if i_stimtype == 2, stimtype = 'NSEM'; end
        inputs.stim_type    = stimtype;
        [StimulusPars DirPars] = Directories_Params_v23(exp_nm, inputs.stim_type, inputs.map_type);
        DirPars.organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', inputs);
        SPars = StimulusPars.slv; 
        if i_stimtype == 1, SPars_WN   = SPars; DirPars_WN   = DirPars; end
        if i_stimtype == 2, SPars_NSEM = SPars; DirPars_NSEM = DirPars; end
    end
    
    %%
    for i_celltype = 1:2    
        if i_celltype == 1; cells = allcells{i_exp}.ONP;  celltype = 'ONPar'; end
        if i_celltype == 2; cells = allcells{i_exp}.OFFP; celltype = 'OFFPar'; end
            %%
        for i_cell = 1:length(cells)
            %%
            clf
            cid = cells(i_cell);
            cell_savename = sprintf('%s_%d', celltype,cid);
            display(sprintf('Working on %s  Cell: %s', exp_nm,cell_savename));
            for i_stimtype  = 1:2
   
                    if i_stimtype == 1, stimtype = 'WN';   SPars = SPars_WN;   DirPars = DirPars_WN;   subplot(7,1,[1 2 3]); hold on; end
                    if i_stimtype == 2, stimtype = 'NSEM'; SPars = SPars_NSEM; DirPars = DirPars_NSEM; subplot(7,1,[5 6 7]); hold on; end
                    set(gca,'fontsize',10)
                    params.evalblocks = SPars.TestBlocks;
                    params.bins = params.bins_per_frame *length(SPars.testframes); 
                    eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', DirPars.organizedspikesdir, cell_savename));

                    blocks  = length(params.evalblocks);                
                    for i_blk = 1 : blocks
                        blknum = params.evalblocks(i_blk);
                        spikes = organizedspikes.block.t_sp_withinblock{blknum};
                        spikes = spikes(find(spikes <=plotsecs));
                        spikes = spikes(find(spikes > 0 ) );
                        yval = i_blk*ones(size(spikes));
                        plot(spikes,yval,'k.','markersize',MS)
                    end
                    metrics = raster_scores{i_exp}.scores.stim_type{i_stimtype}.celltype{i_celltype};
                    spikespersec = metrics.spikespersec(i_cell);
                    trial_variance = metrics.trial_variance(i_cell);
                    temporal_variance = metrics.temporal_variance(i_cell);
                    uop_bps = metrics.bits_per_spike(i_cell);

                    xlim([0,plotsecs]); ylim([1,blocks]);
                    xlabel('seconds')
                    set(gca,'xtick',1:plotsecs)
                    
                    title(sprintf('%s %s %d: {Hz = %1.1e, TrialVariance = %1.1e, SignalVariance = %1.1e, UOP BPS = %1.1e}',...
                        stimtype, celltype, cid, spikespersec,trial_variance,temporal_variance,uop_bps));

            end
            orient landscape
            eval(sprintf('print -dpdf %s/rast_met_%s', plotdir, cell_savename))
        end
        
        
    end
end
end   
%}