% calculate PCA-STCs 


% White Noise
% Get time-stim data

WN_datafile = '2015-02-24-5/streamed/data006/data006';
datarun=load_data(WN_datafile)
datarun=load_params(datarun)

datarun.names.rrs_neurons_path='/Volumes/Analysis/2015-02-24-5/streamed/data006/data006.neurons';
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);
datarun=load_data(datarun,opt);

triggers=datarun.triggers; %onsets of the stimulus presentation

mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/BW-20-8-0.48-11111-16x16.xml';

 [mov,height,width,duration,refresh] = get_movie_ath(mdf_file,...
    triggers, 1,2);

[mvi] = load_movie(mdf_file, triggers);

% 'refresh' parameter is the amount of time a particular image
% is shown before it is changed.
subdivideRefresh= 8; % So, each frame refresh is subdivided into this many parts. 


% For example , if you are showing BW-10-4 .. The refresh for the movie
% would be 8.33*4 ms .. But we still want to have temporal resolution ofs
% 8.33ms. So, set subdivideRefresh=4. If you want temporal resolution of
% 4.16ms, then set subdivideRefresh=8. subdivideRefresh must be an integer

distinctImagesbetweenTriggers=(mean(abs(diff(triggers)))*1000/refresh);
% 'distinctImagesbetweenTriggers' is the number of distinct images which
% are shown between refreshes. 

fr=[];
for itrig=1:length(triggers)
fr=[fr,triggers(itrig)*1000+(refresh/subdivideRefresh)*[0:distinctImagesbetweenTriggers*subdivideRefresh-1]];
end
fr=fr';
% fr are times in ms, when a new bin should start. At times when the movie
% frame changes, a new bin should also start. So, for each bin, there is
% only one movie frame. 
idx=0*fr;
idx([1:8:length(fr)])=1;
frMov = fr(logical(idx));

visCID = 3919;
matlab_ids = find(datarun.cell_ids==visCID);
interestingCellIds=matlab_ids;%[1:length(datarun.cell_ids)];% These would be matlabcellIds - NOT.. datarun.cell_ids;% This would change in future, I think
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
spks(ibin)=sum(double(spikes>=fr(ibin) & spikes<fr(ibin+1)));
end
spk_coll{icell}=spks;
end


% Make movie  
mov=zeros(height,width,3,length(fr)-1);
for ibin=1:length(fr)-1
frameNo=find(fr(ibin)>=frMov,1,'last')-1; 
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

mov=mov-0.5;

[STA,STC] = calculate_sta_stc(spk_coll{1},mov,30);




%% Long Original
% Get time-stim data
rawMovFrames=22080/(8);
[stim,height,width,header_size] = get_raw_movie('/Volumes/Data/2015-02-24-5/Visual/pc2015_02_24_5_data006/20.rawMovie',rawMovFrames,1);
subtract_movies{3}=mean(stim,1);
subtract_movies{3}=mean(stim,1)*0+127.5;
movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
movie=movie/255;
movie= permute(movie,[2,3,1]);

datafile = '2015-02-24-5/data011-from-data006_s_nps/data011-from-data006_s_nps';
datarun=load_data(datafile)

datarun.names.rrs_neurons_path='/Volumes/Analysis/2015-02-24-5/data011-from-data006_s_nps/data011-from-data006_s_nps.neurons';
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);
datarun=load_data(datarun,opt);

triggers=datarun.triggers; %onsets of the stimulus presentation

subdivideRefresh= 1; % So, each frame refresh is subdivided into this many parts. 

distinctImagesbetweenTriggers=100;


% frame intervals
nConditions=2;condDuration = 1.5*60 +2; 
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuronPath);
TTL=double(neuronFile.getTTLTimes());
rawMovFrames = nConditions*condDuration*120;
TTLperTrial=floor(rawMovFrames/100)+1;
nTrials=floor(length(TTL)/TTLperTrial);

samplesPerCondition=condDuration*20000;
ConditionStartTimes=[0:nConditions]*samplesPerCondition+1;

fr=[];
frTTL=[];
TTLstartOld=0;
for iTrial=1:nTrials
iTrial
TTLstart=TTL((iTrial-1)*TTLperTrial+1);
TTLend=TTLstart+condDuration*nConditions*20000;%TTL((iTrial)*TTLperTrial+1);
(TTLstart-TTLstartOld)/20000
it=0;
for itrig=triggers(triggers*20000<TTLend & triggers*20000>=TTLstart)'
fr=[fr,itrig*1000+((1000/120)/subdivideRefresh)*[0:100-1]];

it=it+1;
    if(it==1)
        frTTL=[frTTL,1,zeros(1,99)];
    else
        frTTL=[frTTL,zeros(1,100)];
    end
end
TTLstartOld=TTLstart;

end
fr=fr';
frTTL=frTTL';
frMov=fr;



% Make movie  
% Wrong ?
condition_log=[];
frameNo_log=[];
mov=zeros(height,width,3,length(fr)-1);
idx=1:length(frTTL);idx=idx';
for ibin=1:length(fr)-1
    
    lastframe = find(frTTL==1 & idx<=ibin ,1,'last');
frameNo=floor((ibin-lastframe )/8) +1 ; 
if(frameNo<=size(movie,3))

 mov(:,:,1,ibin)= movie(:,:,frameNo);
 mov(:,:,2,ibin)= movie(:,:,frameNo);
 mov(:,:,3,ibin)= movie(:,:,frameNo);
frameNo_log=[frameNo_log;frameNo];

if(frameNo<size(movie,3)/2)
condition_log=[condition_log;1];
else
condition_log=[condition_log;0];
end

else

 mov(:,:,1,ibin)= 0*movie(:,:,1);
 mov(:,:,2,ibin)= 0*movie(:,:,1);
 mov(:,:,3,ibin)= 0*movie(:,:,1);
 frameNo_log=[frameNo_log;-1];
 
condition_log=[condition_log;-1];

end
end



visCID = 3919;
matlab_ids = find(datarun.cell_ids==visCID);
interestingCellIds=matlab_ids;%[1:length(datarun.cell_ids)];% These would be matlabcellIds - NOT.. datarun.cell_ids;% This would change in future, I think
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
spks(ibin)=sum(double(spikes>=fr(ibin) & spikes<fr(ibin+1)));
end
spk_coll{icell}=spks;
end


mov=mov(:,:,:,1:length(spk_coll{1}));
condition_log=condition_log(1:length(spk_coll{1}));
[STA0,STC0] = calculate_sta_stc(spk_coll{1}(condition_log==0),mov(:,:,:,condition_log==0),30);
[STA1,STC1] = calculate_sta_stc(spk_coll{1}(condition_log==1),mov(:,:,:,condition_log==1),30);