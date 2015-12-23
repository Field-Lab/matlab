
startup_null_analyse_tenessee

% Analysis of white noise - reSTA
% Dataset 2014-09-10-2/data004

cellList=[6392,621,2603,1758,2582,2375];
targetCell=6393;% 6393

%% Load cells
datafile = '2014-09-10-2/data004';


datarun=load_data(datafile)
datarun=load_sta(datarun)
datarun=load_params(datarun)
datarun=load_ei(datarun,'all','array_type',519);
%get_cell_ids(datarun,type_name) % Vision IDs - used for loading fitted
%STAs!
n_cell=length(datarun.cell_ids);


%% Load STAS
stas=datarun.stas.stas;

vision_id=targetCell;
idx=[1:length(datarun.cell_ids)];
matlab_id=idx(datarun.cell_ids==vision_id);
stas=stas{matlab_id};



    st_temp=zeros(size(stas,2),size(stas,1),1,size(stas,4)); % DOUBT .. Could be a reason for things to fail!!!!!
    for itime=1:30
        st_temp(:,:,:,itime)=mean(stas(:,:,:,end-itime+1),3)'; % DOUBT .. Could be a reason for things to fail!!!!!
    end
    stas_new=st_temp;

stas=stas_new;


%% 
vision_id=targetCell;
idx=[1:length(datarun.cell_ids)];
matlab_id=idx(datarun.cell_ids==vision_id);
cell_ei=datarun.ei.eis{matlab_id};

cell_ei_mag=sum(cell_ei.^2,2);
[v,center_elec]=max(cell_ei_mag);

% Get all cells which are strong at calculated center elctrode
% otherCellsMatlabID;
magCenterElec=[];
for icell=1:length(datarun.ei.eis)
    x=datarun.ei.eis{icell}(center_elec,:);
    magCenterElec=[magCenterElec;sum(x.^2)];
end
figure;stem(magCenterElec);

% get surround electrodes 
centerElecPos=datarun.ei.position(center_elec,:);
distCenterElec = datarun.ei.position-repmat(centerElecPos,[size(datarun.ei.position,1),1]);
distCenterElec=sum(distCenterElec.^2,2);
[sortedDist,nearestChannels]=sort(distCenterElec,'ascend');
elecList=nearestChannels(1:7);

% Get Initial waveform estimates
% For center electrode
init_waveform=datarun.ei.eis{matlab_id}(center_elec,:);
figure;plot(init_waveform);

init_waveform_all_channels=datarun.ei.eis{matlab_id}(elecList,:);
figure;
plot(init_waveform_all_channels');

%% Get raw data from data 
% data014

%getTrialData('/Volumes/Data/2014-08-20-0/data003',center_elec,330,30)
samplingRate=20000;
raw_data_path = '/Volumes/Data/2014-09-10-2/data014';
rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(raw_data_path);
%time = 330; %seconds
display('Analysing Data');
time=15*60;

% get Rawdata for center electrode

data=zeros(time*samplingRate,1);
triggerData=0*data;
allElecData=zeros(time*samplingRate,length(elecList));
cnt=1;   
for j=0:1:time-1 % Depends on size of file
        j
        m1=rawFile.getData(j*samplingRate, samplingRate);

       data(cnt:cnt+1*samplingRate-1,:)=   m1(:,center_elec+1); 
       triggerData(cnt:cnt+1*samplingRate-1,:)=m1(:,0+1);
       allElecData(cnt:cnt+1*samplingRate-1,:)=m1(:,elecList+1); 
       cnt=cnt+20000;
end

    data=data';
    data=double(data);
    
    triggerData=triggerData';
    triggerData=double(triggerData);
    
    allElecData=allElecData';
    allElecData=double(allElecData);

rawFile.close();


dataCenter=data;

%% 
spike_sort_CBP

  %% Load Raw movie ? 
  
[stim,height,width,header_size] = get_raw_movie('/Volumes/Data/2014-09-10-2/Visual/13.rawMovie',120*(15*60+2),1);
%subtract_movies{6}=mean(stim,1);
subtract_movies{1}=mean(stim,1)*0+127.5;
movies{1}=stim-repmat(subtract_movies{1},[120*(15*60+2),1,1]);


cl=[min(movies{1}(:)),max(movies{1}(:))];
figure;
for itime=1:120*12
imagesc(squeeze(movies{1}(itime,:,:)));
colormap gray
caxis([-127.5,127.5]);
colorbar
pause(1/120);
end
%%
%Load spike sorted data
load('/Volumes/Analysis/nishal/cell6393_3.mat');
 
figure;
  imov=1;
  n_cell=6
  for icell=1:n_cell
      subplot(4,2,icell);
      plot(waveforms{imov}{icell});
      title(sprintf('Cell : %d',icell));
  end
  
  subplot(4,2,8);
  plot(init_waveform_all_channels');
  title('Vision Waveform');
  
  correct_cell=1;
  display('Correct cell is 1 ,but it has bimodal amplitude histogram');
  figure;
  hist(spike_amps{1}{correct_cell},20)
  title('Spike amplitude histogram as found from CBP')
  
  display('Will get spikes after Thresholding')
  spike_th=0.2;
  
  spk_tms=(spike_times{imov}{correct_cell}(spike_amps{imov}{correct_cell}>spike_th)); % in samples
  
  

%
idx=[1:numel(triggerData)];
triggerSamples=idx(triggerData<-0.01);

trigSamples=[];
trigSamples(1)=triggerSamples(1);
lastTrig=trigSamples(1);
for iLen=2:numel(triggerSamples)
    if(abs(triggerSamples(iLen)-lastTrig)>20);
    trigSamples=[trigSamples;triggerSamples(iLen)];
    lastTrig=trigSamples(end);
    end
end

clear lastTrig triggerSamples


staReEst=cell(1,1);

% calculate new STA
sampleToFrame=120/20000;
binnedSpks=zeros((15*60+2)*120,1);
staLen=30;
for imov=1
width=size(movies{imov},2);
height=size(movies{imov},3);
staReEst{imov} = zeros(staLen,width,height);
spkCount=0;

    startSample = round((staLen+1)*20000/120);
    for ispk=1:numel(spk_tms)
            spkCount=spkCount+1;
            spkSample=(spk_tms(ispk))
            numTriggers=sum(double(trigSamples<spkSample));
            lastTrig=max(trigSamples(trigSamples<spkSample));
            startFrameNo=floor((numTriggers-1)*100+(spkSample-lastTrig)*sampleToFrame+1)-1 % DOUBT
            
            binnedSpks(startFrameNo)=binnedSpks(startFrameNo)+1;
            if(startFrameNo>staLen)
            for i=1:staLen
            staReEst{imov}(i,:,:) = staReEst{imov}(i,:,:)+movies{imov}(startFrameNo-i+1,:,:);
            end
            end
            
    end
 end

staReEst{imov}=staReEst{imov}/spkCount;

figure;
for itime=1:staLen
    itime
    subplot(1,2,1);
imagesc(squeeze(staReEst{imov}(itime,:,:)));
colormap gray
colorbar
caxis([min(staReEst{imov}(:)),max(staReEst{imov}(:))]);
axis image
%caxis(cl);

 subplot(1,2,2);
imagesc(stas(:,:,itime));
colormap gray
colorbar
caxis([min(stas(:)),max(stas(:))]);
%caxis(cl);
axis image
pause
end


 %% STC?
 addpath('../code_stc/');
 if(vision_id==6393)
 xId = [5:20];
 yID=[5:20];
 end
 
 % 6572
 xID=[10:20];
 yID=[10:20];
 
 mov_stc=zeros(size(movies{1},1),numel(xID)*numel(yID));
 for itime=1:size(movies{1},1)
 xx=squeeze(movies{1}(itime,:,:));
xx=xx(xID,yID);

 mov_stc(itime,:)=xx(:)';
 end
 [STA,STC] = simpleSTC(mov_stc, binnedSpks, 30);
 STAr=reshape(STA,[30,numel(xID),numel(yID)]);
 %[u,s,v] = svds(STC,20);  
 

 figure;
for itime=1:staLen
    itime
    subplot(1,2,1);
imagesc(squeeze(STAr(itime,:,:)));
colormap gray
colorbar
caxis([min(STAr(:)),max(STAr(:))]);
axis image
%caxis(cl);

 subplot(1,2,2);
imagesc(stas(:,:,itime));
colormap gray
colorbar
caxis([min(stas(:)),max(stas(:))]);
%caxis(cl);
axis image
pause
end

 
 %%
%  neuronPath ='/Volumes/Analysis/2014-09-10-2/data014/data014.neurons';
% neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuronPath);
% 
% spks=neuronFile.getSpikeTimes(543) 542/543 in data014 and data015 is 622
% and 621  in data 004 new and old respectively
% spks=neuronFile.getSpikeTimes(6443) 6443 in data014 and data015 is 6393
% and 6392 data 004 new and old respectively

%%
neuronPath ='/Volumes/Analysis/2014-09-10-2/data014/data014.neurons';
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuronPath);
neuron_id_list = neuronFile.getIDList();
TTLTimes=neuronFile.getTTLTimes();
vision_spks=neuronFile.getSpikeTimes(2633); % 6443 Vision spks use neuron ID: 6572 to show that STA code works!
vision_spks=double(vision_spks); 
spk_tms=vision_spks(vision_spks<(15*60)*20000);

 
%
idx=[1:numel(triggerData)];
triggerSamples=idx(triggerData<-0.01);

trigSamples=[];
trigSamples(1)=triggerSamples(1);
lastTrig=trigSamples(1);
for iLen=2:numel(triggerSamples)
    if(abs(triggerSamples(iLen)-lastTrig)>20);
    trigSamples=[trigSamples;triggerSamples(iLen)];
    lastTrig=trigSamples(end);
    end
end

clear lastTrig triggerSamples


staReEst=cell(1,1);

% calculate new STA
sampleToFrame=120/20000;
binnedSpks=zeros((15*60+2)*120,1);
staLen=30;
for imov=1
width=size(movies{imov},2);
height=size(movies{imov},3);
staReEst{imov} = zeros(staLen,width,height);
spkCount=0;

    startSample = round((staLen+1)*20000/120);
    for ispk=1:numel(spk_tms)
            spkCount=spkCount+1;
            spkSample=(spk_tms(ispk))
            numTriggers=sum(double(trigSamples<spkSample));
            lastTrig=max(trigSamples(trigSamples<spkSample));
            startFrameNo=round((numTriggers-1)*100+(spkSample-lastTrig)*sampleToFrame+1) % DOUBT
            
            binnedSpks(startFrameNo)=binnedSpks(startFrameNo)+1;
            if(startFrameNo>staLen)
            for i=1:staLen
            staReEst{imov}(i,:,:) = staReEst{imov}(i,:,:)+movies{imov}(startFrameNo-i+1,:,:);
            end
            end
            
    end
 end

staReEst{imov}=staReEst{imov}/spkCount;

figure;
for itime=1:staLen
    itime
    subplot(1,2,1);
imagesc(squeeze(staReEst{imov}(itime,:,:)));
colormap gray
colorbar
caxis([min(staReEst{imov}(:)),max(staReEst{imov}(:))]);
axis image
%caxis(cl);

 subplot(1,2,2);
imagesc(stas(:,:,itime));
colormap gray
colorbar
caxis([min(stas(:)),max(stas(:))]);
%caxis(cl);
axis image
pause
end

