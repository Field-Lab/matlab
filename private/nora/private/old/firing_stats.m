function [spikes] = firing_stats(exps,stimtypes,celltypes)

window_size = 1000;
time_bins = 30*1000; %30 seconds times 1000 ms
cutoff = 0.5;
n_blocks = 59;


% Load core directories and all eligible cells
BD = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection));

for i_exp = exps    
    for i_stimtype = stimtypes
        % Load master datarun, bookkeep
        exp_nm  = allcells{i_exp}.exp_nm;
        expname = allcells{i_exp}.expname;
        eval(sprintf('load %s/%s/datarun_master.mat', BD.BlockedSpikes,exp_nm));
        if i_stimtype == 1, stimtype = 'WN';   end
        if i_stimtype == 2, stimtype = 'NSEM'; end
        [StimulusPars, exp_info] = StimulusParams(exp_nm, stimtype, 'mapPRJ');  
        secondDir.exp_nm    = exp_nm; 
        secondDir.map_type  = 'mapPRJ';
        secondDir.stim_type = stimtype;
        secondDir.fitname   = 'whatever';
        Dirs.organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', secondDir);
        for i_celltype = celltypes
            % Choose which subset of cells to run
            if i_celltype == 1; cellgroup = allcells{i_exp}.ONP;  celltype = 'ONPar'; end
            if i_celltype == 2; cellgroup = allcells{i_exp}.OFFP; celltype = 'OFFPar'; end
            candidate_cells = [allcells{i_exp}.ONP allcells{i_exp}.OFFP];
            cellgroup = intersect(candidate_cells, cellgroup);
            spikes = zeros(length(cellgroup),n_blocks,time_bins);
            for i_cell = 1:length(cellgroup)
                cid = cellgroup(i_cell);
                cell_savename = sprintf('%s_%d', celltype,cid);
                % Load Blocked-Spikes from preprocessing
                eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', Dirs.organizedspikesdir, cell_savename));
                % Process spikes for glm_execute with proper subroutines
                fitspikes = organizedspikes.block.t_sp_withinblock(2:2:end);
                for i_block = 1:n_blocks
                    spikes_ms = round(1000*fitspikes{i_block});
                    spikes(i_cell, i_block, spikes_ms) = 1;
                end
            end
            slider = zeros(1,1,window_size);
            slider(1,1,:) = gausswin(window_size);
            spike_rate = convn(spikes, slider, 'same');
        end
    end
end
end

function spikesconcat      = subR_concat_fitspikes_fromorganizedspikes(blockedspikes, FitPars)
% AKHeitman 2014-04-14
% Concatenate Spikes from different blocks to a single spike train
% blocekdspikes: needs
%   .t_sp_withinblock
%
% FitPars needs
%   .fittest_skipseconds
%   .tstim
%   .fitframes
%   .FitBlocks


t_start   = 0;
tstim     = FitPars.computedtstim;
FitBlocks = FitPars.FitBlocks;


T_SP = []; blk_count = 0;
dur = tstim * 3600;
for k = FitBlocks
	blk_count = blk_count + 1;
	t_sp_full = blockedspikes.t_sp_withinblock{k} ; % unit of time: sec, 0 for the onset of the block
	t_sp      = t_sp_full(find(t_sp_full >  t_start));
	t_sp = t_sp - t_start;
	t_spcontext = t_sp + ( blk_count -1 )*dur;
	T_SP = [T_SP ; t_spcontext];
end
spikesconcat = T_SP;
end

function raster_spiketimes = subR_createraster(blockedspikes, TestPars)
% AKHeitman 2014-04-14
% Make a raster which takes into account GLM processing
% blocekdspikes: needs
%   .t_sp_withinblock
%
% TestPars needs
%   .fittest_skipseconds
%   .TestBlocks

rasterblocks = TestPars.TestBlocks;
t_start      = TestPars.fittest_skipseconds;

raster_spiketimes = cell(length(rasterblocks),1);

for i_blk = 1 : length(rasterblocks)
	blknum = rasterblocks(i_blk);
	sptimes = blockedspikes.t_sp_withinblock{blknum} - t_start;
	sptimes = sptimes(find(sptimes > 0 ) );
    % HACK NEEDED FOR 2013-10-10-0 and other long runs
    if isfield(TestPars, 'test_skipENDseconds')
        sptimes = sptimes(find(sptimes < (TestPars.test_skipENDseconds - TestPars.fittest_skipseconds - .1)));
    end
    
    raster_spiketimes{i_blk} = sptimes;
end 

end