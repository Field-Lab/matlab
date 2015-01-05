% take a spike train from the same experiment and reformat it

% set these three things
cell='OFFPar_346';
exp='2013-10-10-0';
type='NSEM';

% load CBP spikes here
load('/Volumes/Analysis/nora/CBPSpikeSort/2013-10-10-0/NSEM_dataset_2013-10-10-0_cell346_data001_share.mat')
CBP_spike_train=double(spk_tm)/20000;

% load original organized spikes
load(['/Volumes/Analysis/nora/NSEM/BlockedSpikes/' exp '/' type '_mapPRJ/organizedspikes_' cell '.mat'])

% replace the raw spike times
orig_t_sp=organizedspikes.t_sp;
organizedspikes.t_sp=CBP_spike_train;

% replace the block
orig_block=organizedspikes.block;
for i_block=1:length(orig_block.t_sp)
    clear t_b t_e
    t_b = orig_block.t_frame{i_block}(1);
    t_e = orig_block.t_frame{i_block}(end);
    organizedspikes.block.t_sp{i_block} = organizedspikes.t_sp( (organizedspikes.t_sp > t_b) & (organizedspikes.t_sp < t_e) );
    organizedspikes.block.t_sp_withinblock{i_block} = organizedspikes.block.t_sp{i_block} - t_b;
end

hold on
for i_block=1:5:length(orig_block.t_sp)
    % quick plotting sanity check
    scatter(orig_block.t_sp_withinblock{i_block}, -1+i_block*ones(length(orig_block.t_sp_withinblock{i_block}),1),'.b')
    scatter(organizedspikes.block.t_sp_withinblock{i_block}, 1+i_block*ones(length(organizedspikes.block.t_sp_withinblock{i_block}),1),'.r')
end

save(['/Volumes/Analysis/nora/NSEM/CBPBlockedSpikes/' exp '/' type '_mapPRJ/organizedspikes_' cell '.mat'],'organizedspikes')

