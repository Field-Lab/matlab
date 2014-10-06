
startup_null_analyse_tenessee

% Analysis of white noise - reSTA
% Dataset 2014-09-10-2/data004

cellList=[6392,621,2603,1758,2582,2375];
targetCell=2582;% 6393,622,2601 (bad, very low spike rate, electrode 174), ,2582

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
