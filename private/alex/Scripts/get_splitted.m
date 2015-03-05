function get_splitted(ch,file_list,date,length_expected,prefix,inds)
% pathway=['S:\user\alexandra\Mea\',date,'\MEA_bin\channels\channel',int2str(ch),'\'];
pathway=['F:\',date,'\MEA_bin\channels\channel',int2str(ch),'\'];
% HIGH PASS - for SPIKES
passH=500;%frequency of the high-pass filter
sf=25000;%sampling frequency
fNormH = passH/(sf/2);
orderH=10;%filter order
[aH,bH] = butter(orderH, fNormH, 'high');

% LOW PASS - for LFP
sf=25000;%sampling frequency
passL=180;orderL=2;
fNormL = passL/(sf/2);
[aL,bL] = butter(orderL, fNormL, 'low');

maxCells=2;

NoF=length(file_list);
extrems=zeros(3,0.1*maxCells*length_expected*NoF);
thrs=zeros(NoF,1);
spikes=cell(1,NoF);
lfp=zeros(length_expected,NoF);
last=0;
files=dir(pathway);
for j=1:NoF
    load([pathway, files(file_list(j)+2).name]);
    %getting spikes - each cell for stimulus file; only time stamps
    %of spikes
    dataHPF=filtfilt(aH,bH,data);
    recDuration=length(dataHPF)/25000; % duration of recording in s
    maxSpikesAssumed=100*maxCells*recDuration; % max spike amount to consider    
    
%     tmp=-findpeaks(-dataHPF,'MINPEAKHEIGHT',0,'SORTSTR','descend');
      tmp=[0 diff(sign(diff(dataHPF)))];
    tmp=sort(dataHPF(tmp==2));
    thr=tmp(ceil(maxSpikesAssumed));
%     tmp=dataHPF;
%     tmp(tmp>thr)=0;
%     tmp(tmp<=thr)=-1;
%     tmp=diff(tmp);
tmp=diff(dataHPF>thr);
    timeStamps=find(tmp==-1)+1;

    % remove spikes very close to the beginning and the end
    timeStamps(timeStamps<20)=[];
    timeStamps(timeStamps>length(dataHPF)-30)=[];
    
    thrs(j)=thr;
    waveform=zeros(24,length(timeStamps));
    for k=-10:13
        waveform(k+11,:)=dataHPF(timeStamps+k);
    end
    [extrems(1,last+1:last+length(timeStamps)),ind]=min(waveform(11:end,:));
    ind=ind+10;
    [extrems(2,last+1:last+length(timeStamps)),ind1]=max(waveform(1:11,:));
    extrems(3,last+1:last+length(timeStamps))=ind-ind1; 
    last=last+length(timeStamps);
    spikes{j}=timeStamps/25; % spike times in ms
    clear dataHPF tmp timeStamps waveform
    
    %getting lfp. rows: stimuli; columns: data (each tenth)
%     data = filtfilt(aL, bL, data');
%     data=data(1:25:end);
%     while length(data)<length_expected
%         data=[data; zeros(5000,1)];
%     end
%     lfp(:,j)=data(1:length_expected);
%     clear data
%     if mod(j,100)==0
%         display(j)
%     end
end
extrems=extrems(:,1:last);

save(['F:\',date,'\processing\stat\',prefix,'_',inds,'_CH',int2str(ch),'_',date(5:8)],'extrems','thrs');
name2save=['F:\',date,'\processing\preprocessed\lfp\',prefix,'_',inds,'_CH',int2str(ch),'_',date(5:8),'.mat'];
save(name2save,'lfp');
save(['F:\',date,'\processing\preprocessed\spikes\',prefix,'_',inds,'_CH',int2str(ch),'_',date(5:8),'.mat'],'spikes');
clear lfp spikesBin