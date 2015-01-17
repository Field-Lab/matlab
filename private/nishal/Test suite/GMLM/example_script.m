%%
clear datarun
datarun.names.rrs_neurons_path='/Volumes/Analysis/2012-08-09-3/data002/data002.neurons';
    
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);
datarun=load_data(datarun,opt);

triggers=datarun.triggers; %onsets of the stimulus presentation

%%
mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-8-1-0.48-11111.xml';

 [mov,height,width,duration,refresh] = get_movie_ath(mdf_file,...
    triggers, 1,2);

[mvi] = load_movie(mdf_file, triggers);

% 'refresh' parameter is the amount of time a particular image
% is shown before it is changed.
subdivideRefresh= 10; % So, each frame refresh is subdivided into this many parts. 


% For example , if you are showing BW-10-4 .. The refresh for the movie
% would be 8.33*4 ms .. But we still want to have temporal resolution of
% 8.33ms. So, set subdivideRefresh=4. If you want temporal resolution of
% 4.16ms, then set subdivideRefresh=8. subdivideRefresh must be an integer

distinctImagesbetweenTriggers=round(mean(abs(diff(triggers)))*1000/refresh);
% 'distinctImagesbetweenTriggers' is the number of distinct images which
% are shown between triggers. 

fr=[];
for itrig=1:length(triggers)
fr=[fr,triggers(itrig)*1000+(refresh/subdivideRefresh)*[0:distinctImagesbetweenTriggers*subdivideRefresh]];
end
fr=fr';
% fr are times in ms, when a new bin should start. At times when the movie
% frame changes, a new bin should also start. So, for each bin, there is
% only one movie frame. 

frMov=[];
for itrig=1:length(triggers)
frMov=[frMov,triggers(itrig)*1000+(refresh)*[0:distinctImagesbetweenTriggers]];
end
frMov=frMov';

% frMov are times when the frames change. frMov should times
% should be every subdivideRefresh-th of fr times.



interestingCellIds=[7];%[1:length(datarun.cell_ids)];% These would be matlabcellIds - NOT.. datarun.cell_ids;% This would change in future, I think
nCells=length(interestingCellIds);

spk_coll=cell(nCells,1);
for icell=1:nCells
    icell
    % Next, count spikes in each of these bins.
cellID=interestingCellIds(icell);
spikes=datarun.spikes{cellID};
% spikes are in s. convert to ms
spikes=round(spikes*1000);
spks=zeros(length(fr)-1,1);
for ibin=1:length(fr)-1
    if(mod(ibin,1000)==1)
    ibin
    end
spks(ibin)=sum(double(spikes>=fr(ibin) & spikes<fr(ibin+1)));
end
spk_coll{icell}=spks;
end


% Make movie 
fr2=frMov;
mov=zeros(height,width,3,length(fr2)-1);
for ibin=1:length(fr2)-1
    if(mod(ibin,1000)==1)
    ibin
    end
frameNo=find(fr2(ibin)>=frMov,1,'last')-1; 
if(frameNo<duration)
 F = round(mvi.getFrame(frameNo).getBuffer);
 mov(:,:,1,ibin)= reshape(F(1:3:end),width,height)';
 mov(:,:,2,ibin)= reshape(F(2:3:end),width,height)';
 mov(:,:,3,ibin)= reshape(F(3:3:end),width,height)';
else
    break
end
end

% Correct for short movie
% If there was an error because movie was small in length, correct for it.
finalLen=ibin;
mov=mov(:,:,:,1:finalLen);
for icell=1:nCells
   spk_coll{icell}=spk_coll{icell}(1:finalLen); 
end

% Play movie ? 
figure
for ibin=1:10
imagesc(mov(:,:,:,ibin));
title(sprintf('%d',ibin));
pause(1/120)
end

% have spk_coll and mov at end! 
%% 
trainGMLM(spk_coll{1},mov,2)