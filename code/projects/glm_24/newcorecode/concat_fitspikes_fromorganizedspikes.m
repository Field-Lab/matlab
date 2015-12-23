% AKHeitman 2014-04-14
% StimulusPars.slv  dependent!
% Where we group spikes to fit the model
% Very careful with seconds timing  (use computedtstim!!)


% MMOTE: Make a new version in "true time" 2014-05-03
function spikesconcat = concat_fitspikes_fromorganizedspikes(blockedspikes, FitPars);

% blocekdspikes: needs
%   .t_sp_withinblock

% FitPars needs
%   .fittest_skipseconds
%   .tstim
%   .fitframes
%   .FitBlocks


t_start   = FitPars.fittest_skipseconds;
tstim     = FitPars.computedtstim;
fitframes = FitPars.fitframes;
FitBlocks = FitPars.FitBlocks;


T_SP = []; blk_count = 0;
dur = tstim * length(fitframes);
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
        