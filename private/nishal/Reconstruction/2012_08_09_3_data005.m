
startup_null_analyse_tenessee



datafile = '2012-08-09-3/data002';


opt=struct('load_all',true);
datarun=load_data(datafile,opt)
datarun=load_sta(datarun)
datarun=load_params(datarun)
datarun=load_ei(datarun,'all','array_type',519);
CellID002=datarun.cell_types{1}.cell_ids;

neuronPairsRefVsNew = load('/Volumes/Analysis/2012-08-09-3/data005/cellmatch_mapEI_from_data002.txt');%crossIdentifyNeuronIDs('/Volumes/Analysis/2012-08-09-3/data002', '/Volumes/Analysis/2012-08-09-3/data005',OnParasol_vis_id);
   
CellID005=[];
secCol = neuronPairsRefVsNew(:,2);
for icell=1:length(CellID002)
 CellID005=[CellID005;secCol(find(neuronPairsRefVsNew(:,1)==CellID002(icell)))];
end



%%

clear datarun

datafile = '2012-08-09-3/data005';


opt=struct('load_all',true);
datarun=load_data(datafile,opt)
datarun=load_sta(datarun)
datarun=load_params(datarun)
datarun=load_ei(datarun,'all','array_type',519);
CellID002=datarun.cell_types{1}.cell_ids;

triggers=datarun.triggers; %onsets of the stimulus presentation
refresh=1000/120;%datarun.stimulus.refresh_period;

matlab_ids=zeros(length(CellID005),1);

for icell=1:length(CellID005)
matlab_ids(icell)=find(datarun.cell_ids==CellID005(icell));
end


% 'refresh' parameter is the amount of time a particular image
% is shown before it is changed.
subdivideRefresh= 1; % So, each frame refresh is subdivided into this many parts. 


% For example , if you are showing BW-10-4 .. The refresh for the movie
% would be 8.33*4 ms .. But we still want to have temporal resolution of
% 8.33ms. So, set subdivideRefresh=4. If you want temporal resolution of
% 4.16ms, then set subdivideRefresh=8. subdivideRefresh must be an integer

distinctImagesbetweenTriggers=100%round(mean(abs(diff(triggers)))*1000/refresh);
% 'distinctImagesbetweenTriggers' is the number of distinct images which
% are shown between refreshes. 

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


% TODO - nishal!! CHANGE
% 
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

load('/Volumes/Analysis/nishal/recons_2012_08_09_3_data005_2.mat');

%load_raw_move_data005
%load('/Volumes/Analysis/nishal/recons_2012_08_09_3_data005_mov.mat');
%height=size(X,3);
%width=size(X,2);

% Make movie  
duration=3600;
iTrial=1;
lastTrial=0;
nTrials=floor(length(fr)/(3600*3));
% mov for single pixel only
pixX=160;
pixY=80;
load_raw_move_data005

mov=zeros(length(fr)-1,1);
newTrialStart=zeros(nTrials,1);
newTrialStart(1)=1;
binType=zeros(length(fr)-1,1);

for ibin=1:length(fr)-1
    if(mod(ibin,1000)==1)
    ibin
    end
frameNo=find(fr(ibin)>=frMov,1,'last'); 
frameNoTrial=(frameNo-3600*3*(iTrial-1));
if(frameNoTrial>3600*3)
iTrial=iTrial+1;
frameNoTrial=(frameNo-3600*3*(iTrial-1));
newTrialStart(iTrial)=ibin;
end

if(frameNoTrial<=3600) % Testing
 %F = squeeze(X(frameNo,:,:,:));
 frameNoTest=frameNoTrial;
 F=squeeze(X(frameNoTest))';
 mov(ibin)= F;
 %mov(:,:,2,ibin)= F;
 %mov(:,:,3,ibin)= F;
 binType(ibin)=1; % Testing frame
 else % Training
 frameNoTrain= 3600+(iTrial-1)*7200+(frameNoTrial-3600);
 F=squeeze(X(frameNoTrain))';
 mov(ibin)= F;
 %mov(:,:,2,ibin)= F;
 %mov(:,:,3,ibin)= F;
 binType(ibin)=2; % Training frame
end
end



% Correct for short movie
% If there was an error because movie was small in length, correct for it.
finalLen=ibin;
nCells=length(spk_coll);
mov=mov(1:finalLen);
for icell=1:nCells
   spk_coll{icell}=spk_coll{icell}(1:finalLen); 
end

% % Play movie ? 
% figure
% for ibin=1:100
% imagesc(mov(:,:,1,ibin));
% colormap gray
% title(sprintf('%d',ibin));
% pause(datarun.stimulus.refresh_period/1000)
% end

