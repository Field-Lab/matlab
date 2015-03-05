%% get filtered versions of LFP and spike times
clear all
cd('F:\templates_code')
NoF=750;
mainPath='F:\20110913\MEA_bin\channels\channel';
heka=dir('F:\20110913\HEKA\');
file_types=ones(1,NoF);
file_types(1:150:end)=0;
file_types(150:150:end)=0;

length_expected=7700;
file_list=find(file_types(1:NoF)==1)+1354;
date='20110913';
prefix='back';
inds='1355_2105';
for ch=1:60
    tic
    if ch~=15
        ch
        get_splitted(ch,file_list,date,length_expected,prefix,inds)
    end
    toc
end




clear all
cd('F:\templates_code')
NoF=3004;
mainPath='F:\20110913\MEA_bin\channels\channel';
heka=dir('F:\20110913\HEKA\');
file_types=ones(1,NoF);
for i=1:NoF
    if ~isempty(regexp(heka(i+2).name,'60000', 'once'))
        file_types(i)=0;      % spontan. activity recoring   
    elseif ~isempty(regexp(heka(i+2).name,'sh_2', 'once'))
        file_types(i)=-1;  %nd switch
    else
        file_types(i)=1;    % quick stimulation
    end
end
length_expected=7700;
file_list=find(file_types(1:NoF)==1);
file_list=file_list(file_list>2255);
date='20110913';
prefix='quick';
inds='2255_3004';
for ch=1:60
    tic
    if ch~=15
        ch
        get_splitted(ch,file_list,date,length_expected,prefix,inds)
    end
    toc
end











a=1976167 
k=0;
for i=1:length(spikes)
    k=k+length(spikes{i});
    if k>a
        i
        a=10000000000;
    end
end







clear all

% GET READY
cd('F:\20110913\processing')

% get variables for standard plot - channelwise
NoF=740; % number of cells with data (148 per ND)
stim=[2951.2, 4935.3]; % stimulus times in ms
stimBin=floor(stim);
stLength=stim(2)-stim(1);
stBinLength=stimBin(2)-stimBin(1);
stimLFP=floor(stim);
stLFPLength=stimLFP(2)-stimLFP(1);

for i=1:60
    if i~=15
        load(['F:\20110913\processing\preprocessed\lfp\quick_2255_3004_CH',int2str(i),'_0913'])
        load(['F:\20110913\processing\preprocessed\spikes\quick_2255_3004_CH',int2str(i),'_0913'])        
        OnSpikes=zeros(3,NoF); % max (row 1), its std (row 2) and its LP (row 3) of averaged firing rate to the onset of stimulus
        OffSpikes=zeros(3,NoF); % max (row 1), its std (row 2) and its LP (row 3) of averaged firing rate to the offset of stimulus (~1 sec after flash end)
        OnLFP=zeros(3,NoF);
        OffLFPminima=zeros(3,NoF);
        OffLFPmaxima=zeros(3,NoF);
        spontSpikes=zeros(2,NoF);% mean (row 1) and standard deviation (row 2) of spontaneous spikes;
        tic
        for k=1:2:NoF
            conv=convolved(spikes{k},40,7700);
            datplotSpikes=conv(stimBin(1):stimBin(2)+stBinLength);
            
            datplotLFP=lfp(stimLFP(1):stimLFP(2)+stLFPLength,k);
            
            [OnSpikes(1,k),OnSpikes(3,k)]=max(datplotSpikes(1:stimBin(2)-stimBin(1)));
            [OffSpikes(1,k),OffSpikes(3,k)]=max(datplotSpikes(stimBin(2)-stimBin(1):end));
            spontSpikes(1,k)=mean(conv(1:floor(stimBin(1)/2)));
            spontSpikes(2,k)=std(conv(1:floor(stimBin(1)/2)));
            [OnLFP(1,k),OnLFP(3,k)]=min(datplotLFP(1:stimLFP(2)-stimLFP(1)));
            [OffLFPminima(1,k),OffLFPminima(3,k)]=min(datplotLFP(stimLFP(2)-stimLFP(1):end));
            [OffLFPmaxima(1,k),OffLFPmaxima(3,k)]=max(datplotLFP(stimLFP(2)-stimLFP(1):end));
            
            OnSpikes(2,k)=std(datplotSpikes(1:stimBin(2)-stimBin(1)));
            OffSpikes(2,k)=std(datplotSpikes(stimBin(2)-stimBin(1):end));
            OnLFP(2,k)=std(datplotLFP(1:stimLFP(2)-stimLFP(1)));
            OffLFPminima(2,k)=std(datplotLFP(stimLFP(2)-stimLFP(1):end));
            OffLFPmaxima(2,k)=std(datplotLFP(stimLFP(2)-stimLFP(1):end));
        end
        for k=2:2:NoF
            conv=convolved(spikes{k},40,7700);
            datplotSpikes=conv(stimBin(1):stimBin(2)+stBinLength);
            
            datplotLFP=lfp(stimLFP(1):stimLFP(2)+stLFPLength,k);
            
            [OnSpikes(1,k),OnSpikes(3,k)]=max(datplotSpikes(stimBin(2)-stimBin(1):end));
            [OffSpikes(1,k),OffSpikes(3,k)]=max(datplotSpikes(1:stimBin(2)-stimBin(1)));
            spontSpikes(1,k)=mean(conv(1:floor(stimBin(1)/2)));
            spontSpikes(2,k)=std(conv(1:floor(stimBin(1)/2)));
            [OnLFP(1,k),OnLFP(3,k)]=min(datplotLFP(stimLFP(2)-stimLFP(1):end));
            [OffLFPminima(1,k),OffLFPminima(3,k)]=min(datplotLFP(1:stimLFP(2)-stimLFP(1)));
            [OffLFPmaxima(1,k),OffLFPmaxima(3,k)]=max(datplotLFP(1:stimLFP(2)-stimLFP(1)));
            
            OnSpikes(2,k)=std(datplotSpikes(stimBin(2)-stimBin(1):end));
            OffSpikes(2,k)=std(datplotSpikes(1:stimBin(2)-stimBin(1)));
            OnLFP(2,k)=std(datplotLFP(stimLFP(2)-stimLFP(1):end));
            OffLFPminima(2,k)=std(datplotLFP(1:stimLFP(2)-stimLFP(1)));
            OffLFPmaxima(2,k)=std(datplotLFP(1:stimLFP(2)-stimLFP(1)));
        end
        toc
        save(['F:\20110913\processing\ready\ready_quick_2255_3004_CH',int2str(i),'_0913.mat'],'OnSpikes','OffSpikes','OnLFP','OffLFPminima','OffLFPmaxima','spontSpikes');
        clear spikes lfp
    end
end
