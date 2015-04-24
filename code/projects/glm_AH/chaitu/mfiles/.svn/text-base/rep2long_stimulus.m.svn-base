function [sptimes_aligned Xlong] = rep2long_stimulus(Xrep,ntrials,trial_duration,nrepframes,sptimes,dt,blanklevel)
% Arguments:

% Xrep - repeat stimulus
% ntrials - number of trials
% trial_duration - trial duration
% nrepframes  - number of frames per repeat trail
% sptimes - list of analog spike_times for whole repeat experiment
% dt - stimulus frame time step
% blanklevel - level at which to set blank frames
% Returns:

% Xlong - stimulus containing repeat trials with blank buffer frames
% sptimes_aligned - realigned spike times


buffertime = trial_duration - nrepframes*dt;
nblanks = floor(buffertime/dt);

nrepframes_eff = nrepframes + nblanks;
trial_duration_eff = nrepframes_eff*dt;

slack = trial_duration - trial_duration_eff;

if (nargout > 1)
    Xlong = zeros(size(Xrep,1),nrepframes_eff*ntrials);
end

sptimes_aligned = [];


fprintf('Frame duration is %0.5f sec (slack is %0.5f). Each trial has %d frames (%d stimulus, %d blanks)\n',trial_duration_eff,slack,nrepframes_eff,nrepframes,nblanks);


for j=1:ntrials
    
    if (nargout > 1)
        trial_idx = (j-1)*nrepframes_eff+1:(j-1)*nrepframes_eff+nrepframes;
        blank_idx = (j-1)*nrepframes_eff+nrepframes+1:(j-1)*nrepframes_eff+nrepframes+nblanks;
        Xlong(:,trial_idx) = Xrep;
        Xlong(:,blank_idx) = ones(size(Xrep,1),nblanks).*blanklevel;
    end
    
    % Get all spikes in this trial period
    starttime = (j-1)*trial_duration;
    finishtime = starttime+trial_duration_eff;
    sp_idx = sptimes > starttime & sptimes <= finishtime;
    sptimes_aligned = [sptimes_aligned; (sptimes(sp_idx)-slack*(j-1))];
end



    