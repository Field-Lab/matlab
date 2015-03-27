% AKHeitman 2014-11-26
% Hack Code is done and already computed with
% sigmas = [.5 1 2 4 7 10 15 20 25 35 50 70 100 150 250 500 1000];

% Main Call: rasterprecision_fullxval
% Computes BPS of raster given average firing rate .. full XVAL
% Computes Fraction Error given average firing rate .. full XVAL

% Saved into raster_metrics_temporalprecision folder within BlockedSpikes

%%

clear; close all; clc

% SET DIRECTORIES
BD = NSEM_BaseDirectories;
DirPars.rastmetdir = sprintf('%s/raster_performance', BD.Cell_Selection);
if ~exist(rastmetdir, 'dir'), mkdir(rastmetdir); end
DirPars.rast_dir = BD.BlockedSpikes;

% LOAD CELL NUMBERS / INITIALIZE SOLUTION STRUCTURE
raster_scores = allcells;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection));  % allcells structire



% TIME SCALES
sigmas = [.5 1 2 4 7 10 15 20 25 35 50 70 100 150 250 500 1000];
params.bindur     = .00083275;
params.bins_per_frame = 10;
metricparams.smoothbins    = sigmas;


% WHICH CELLS TO RUN THROUGH
exps = 1:4;
stimtypes = 1:2;
celltypes = 1:2;
%% LOOP THROUGH AND COMPUTE RASTER SCORES
for i_exp = exps
    for i_stimtype = stimtypes   
        for i_celltype = celltypes
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
            
            d_save = sprintf('%s/raster_metrics_temporalprecision', DirPars.organizedspikesdir);
            
            if ~exist(d_save,'dir'), mkdir(d_save); end
            
          
            %%
            for i_cell = 1:length(cells)
                %% 
                clear raster_precision
                % LOAD CELLS INFORMATION AND GENERATE RASTER %
                cid = cells(i_cell);
                cell_savename = sprintf('%s_%d', celltype,cid);
               % display(sprintf('Working on %s  Cell: %s', exp_nm,cell_savename));
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
                display(sprintf('Working on %s  Cell: %s', exp_nm,cell_savename));
                
                raster_precision = rasterprecision_fullXVAL_old(logicalspike,sigmas);
                raster_precision.cid = cid;
                raster_precision.exp_nm = exp_nm;
                raster_precision.celltype = celltype;
                raster_precision.name  = cell_savename;
                
                eval(sprintf('save %s/temporalprecision_%s.mat raster_precision', d_save,cell_savename));


                
            end     
        end
    end
end
