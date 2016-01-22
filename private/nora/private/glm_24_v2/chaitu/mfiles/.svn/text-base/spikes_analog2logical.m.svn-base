function D = spikes_analog2logical(sp,nrepeatframes,frame_dt,fac,fulltrial_duration,ntrials)

trial_duration = nrepeatframes*frame_dt; % duration of the repeated trial (in sec)
blankframes = 0;%ceil((fulltrial_duration - trial_duration)/frame_dt); % number of blank frames appended after nrepeat frames

%spT = nrepeatframes*fac; % number of timebins at spike resolution for repeat trial
spTeff = (nrepeatframes+blankframes)*fac;%floor(fulltrial_duration/trainpars.dt);
D = logical(false(spTeff,ntrials)); % Matrix of spikes
    
% Bin the spikes for each trial to create a binary matrix (time x trial)
for k2=1:ntrials
    %sp_blocked{k} = sp(sp > (k-1)*5 & sp < k*5);
    spchunk = sp(sp > (k2-1)*fulltrial_duration & sp < k2*fulltrial_duration) - (k2-1)*fulltrial_duration; % get all spikes in this 5sec period
    D(:,k2) = double(binspikes(spchunk,frame_dt/fac, frame_dt/fac*spTeff))'; %only count spikes in [0,nrepeatframes*frame_dt]
end


