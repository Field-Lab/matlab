startup

%%
%plot rawData%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tstart=0;%pick time

tstep=1;

stop=.2;


samplingrate=20000;
%initiate data
raw_data_path='/Volumes/Data/2014-04-15-4/data001';
rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(raw_data_path);

%%
% HIGH PASS - for SPIKES
passH=500;%frequency of the high-pass filter
sf=20000;%sampling frequency
fNormH = passH/(sf/2);
orderH=10;%filter order
[aH,bH] = butter(orderH, fNormH, 'high');
clear passH fNormH orderH

%%
spikes=zeros(24*7,35000);
spikeTimes=zeros(1,350000);
maxSpikesSecond=80;
last=0; % last spike processed
tic
for j=0:60:899 % Depends on size of file
    j
    data=zeros(1200000,7);
    cnt=1;
    for i=j:j+59
        m1=rawFile.getData(i*samplingrate, samplingrate);
        data(cnt:cnt+19999,1)=m1(:,1);
        data(cnt:cnt+19999,2)=m1(:,2);
        data(cnt:cnt+19999,3)=m1(:,3);
        data(cnt:cnt+19999,4)=m1(:,4);
        data(cnt:cnt+19999,5)=m1(:,5);
        data(cnt:cnt+19999,6)=m1(:,6);
        data(cnt:cnt+19999,7)=m1(:,7);
        cnt=cnt+20000;
    end
    
    

    
    data1=filtfilt(aH,bH,data(:,4)');
 
    recDuration=length(data1)/sf; % duration of recording in s
    maxSpikesAssumed=round(maxSpikesSecond*recDuration); % max spike amount to consider
    
    tmp=[0 diff(sign(diff(data1)))];
    tmp=sort(data1(tmp==2));
    thr=tmp(ceil(maxSpikesAssumed));
    tmp=diff(data1>thr);
    timeStamps=find(tmp==-1)+1;
    
    % remove spikes very close to the beginning and the end
    timeStamps(timeStamps<20)=[];
    timeStamps(timeStamps>length(data)-30)=[];
    
    
    waveform=zeros(24*7,length(timeStamps));
    for k=-9:14
        waveform(k+10,:)=data(timeStamps+k,1);
        waveform(k+10+24,:)=data(timeStamps+k,2);
        waveform(k+10+24*2,:)=data(timeStamps+k,3);
        waveform(k+10+24*3,:)=data(timeStamps+k,4);
        waveform(k+10+24*4,:)=data(timeStamps+k,5);
        waveform(k+10+24*5,:)=data(timeStamps+k,6);
        waveform(k+10+24*6,:)=data(timeStamps+k,7);
    end
    
    timeStamps=timeStamps/sf+j; % spike times in s
    spikes(:,last+1:last+length(timeStamps))=waveform; % waveforms
    spikeTimes(last+1:last+length(timeStamps))=timeStamps;
    last=last+length(timeStamps);
end
toc
spikes(:,last+1:end)=[];
spikeTimes(:,last+1:end)=[];
clear data data1 m1 tmp waveform
save('/home/vision/Dropbox/Lab/Development/matlab-standard/private/nishal/rawdata stuff/el1_2.mat','spikes','spikeTimes')