%% get filtered versions of LFP and spike times
clear all
cd('F:\templates_code')
NoF=750;
date='20110928';
mainPath='S:\user\alexandra\Mea\20110928\MEA_bin\channels\channel';
heka=dir('S:\user\alexandra\Mea\20110928\HEKA\');
file_types=ones(1,NoF);
file_types(1:150:end)=0;
file_types(150:150:end)=0;

length_expected=7700;
file_list=find(file_types(1:NoF)==1);
date='20110928';
prefix='quick';
inds='1_750';
used=[15];
while length(used)<60
    for ch=1:60
        tic
        if isempty(find(used==ch, 1))
            if exist([mainPath,int2str(ch)],'dir')
                ch
                get_splitted(ch,file_list,date,length_expected,prefix,inds)
                used=[used ch];
            end
        end
        toc
    end
end



clear all

% GET READY
cd('F:\20110928\processing')

% get variables for standard plot - channelwise
NoF=821; % number of cells with data (148 per ND)
stim=[2951.2, 4935.3]; % stimulus times in ms
stimBin=floor(stim);
stLength=stim(2)-stim(1);
stBinLength=stimBin(2)-stimBin(1);
stimLFP=floor(stim);
stLFPLength=stimLFP(2)-stimLFP(1);

for i=1:60
    if i~=15
        load(['F:\20110928\processing\preprocessed\lfp\quick_1_832_CH',int2str(i),'_0928'])
        load(['F:\20110928\processing\preprocessed\spikes\quick_1_832_CH',int2str(i),'_0928'])        
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
        save(['F:\20110928\processing\ready\ready_quick_1_832_CH',int2str(i),'_0928.mat'],'OnSpikes','OffSpikes','OnLFP','OffLFPminima','OffLFPmaxima','spontSpikes');
        clear spikes lfp
    end
end
