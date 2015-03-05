%% get filtered versions of LFP and spike times
clear all

NoF=750;
mainPath='F:\20110524\MEA_bin\channels\channel';
cd(strcat(mainPath,'5'));
c=dir;
load(strcat(mainPath,'5\',c(3).name));
len=ceil(length(data)/10);


% HIGH PASS - for SPIKES
passH=500;%frequency of the high-pass filter
sf=25000;%sampling frequency
fNormH = passH/(sf/2);
orderH=10;%filter order
[aH,bH] = butter(orderH, fNormH, 'high');

stim=[4968.7, 6952.8];
stLength=stim(2)-stim(1);
len=7900;% length of data in ms
bin=80;%width in ms
shift=20;%shift in ms; =bin to avoid shifting
steps=(len-bin)/shift+1;
stimBin=floor(stim/shift);
stBinLength=stimBin(2)-stimBin(1);
stimLFP=floor(stim*2.5);
stLFPLength=stimLFP(2)-stimLFP(1);
begins=1:bin:len;
ends=shift:bin:len;

for i=4:60
    tic
    display(i)
    if i~=15
         pathName=[mainPath,int2str(i)];
        cd(pathName);
        fileList=dir;
        spikes=cell(1,NoF);
        j=3
        load(fileList(j+2).name);
        maxSpikesAllowed=length(data)/50;
        %getting spikes - each cell for stimulus file; only time stamps
        %of spikes
        dataHPF = filtfilt(aH,bH, data');
        dataHPF=dataHPF.*dataHPF.*dataHPF;
        spikesTotal=maxSpikesAllowed+1;
        coef=1;
        while spikesTotal>maxSpikesAllowed
            test=dataHPF;
            test1=test;
            test1(test>-coef*std(test([1:10000 end-10000:end])))=0;
            test=diff(test1);
            pks=[0; (test==test1(2:end) & test)];
            spikesTotal=nnz(pks);
            coef=coef+0.5;
        end
        pks=pks(1:len*25);
        tic
        tmp=logical(sum(reshape(pks,25,len)));
%         tmp=logical(tmp);
        for bin_cnt=1:length(begins)
            spikesBin(j,bin_cnt)=sum(tmp(begins(bin_cnt):ends(bin_cnt)));
        end
        toc
        timeStamps=find(pks);
        clear pks test test1 dataHPF
        spikes{j}=timeStamps/25;
        
        load('F:\20110524\processing\preprocessed\1-750_CH1_0524.mat')
        spikesBin=zeros(length(spikes),ceil(steps));
        for k=1:length(spikes)
            for m=1:ceil(steps)
                spikesBin(k,m)=sum(logical((spikes{k}>shift*(m-1)).*(spikes{k}<shift*(m-1)+bin)));
            end
        end
        spikesBin=spikesBin';
        spikes=spikesBin;
        save(['F:\20110524\processing\binned\1-750_CH',int2str(i),'_binned_0524.ready.mat'],'spikesBin');
    end
    toc
end