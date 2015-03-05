mea2mat=[0 21 19 15 16 12 10 0
    24 22 20 17 14 11 9 7
    26 25 23 18 13 8 6 5
    29 30 28 27 4 3 1 2
    32 31 33 34 57 58 60 59
    35 36 38 43 48 53 55 56
    37 39 41 44 47 50 52 54
    0 40 42 45 46 49 51 0]';

dates=['20110131'; '20110208'; '20110210'; '20110303'; '20110306'; '20110308'; '20110311'; '20110314'; '20110315'; '20110729'];

%% Filter parameters
% HIGH PASS - for SPIKES
passH=500;%frequency of the high-pass filter
sf=25000;%sampling frequency
fNormH = passH/(sf/2);
orderH=10;%filter order
[aH,bH] = butter(orderH, fNormH, 'high');
len=60000;% length of data in ms
bin=80;%width in ms
shift=20;%shift in ms; =bin to avoid shifting
steps=(len-bin)/shift+1;
begins=1:shift:len-bin+shift;
ends=bin:shift:len;
for i=1:length(begins)
    intrv(1:bin,i)=(begins(i):ends(i));
end
for dat_cnt=1:10
    date=dates(dat_cnt,:);
    if dat_cnt>5
        path=['F:\',date,'\MEA_bin\channels\channel'];
    else
        path=['F:\',date,'\channels\channel'];
    end
    cd([path,'1']);
    a=dir;
    cnt=1;
    m=[];
    for i=3:length(a)
        if ~isempty(regexp(a(i).name,'spontan'))
            m(cnt)=i;
            cnt=cnt+1;
        end
    end
    spikes=[];
    clear spikesBin
    for i=1:60
        i
        for j=1:length(m)        
            a=dir([path,int2str(i)]);
            load ([path,int2str(i),'\',a(m(j)).name]);
            data=data(1:1500000);
            maxSpikesAllowed=length(data)/50;
            if size(data,1)==1
                data=data';
            end
            dataHPF = filtfilt(aH,bH, data);
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
            % binning: assumes no more than 1 spike per 1 ms. Not for spike times!
            tmp=logical(sum(reshape(pks,25,len)));
            tmp=tmp(intrv);
            spikesBin(i,1:length(begins),j)=sum(tmp,1);
            % get spike times
            timeStamps=find(pks);
            clear pks test test1 dataHPF
            spikes{i,j}=timeStamps/25;
            cnt=cnt+1;
        end
    end
    save(['F:\',date,'\processing\preprocessed\spontan_all_spikes'],'spikes','spikesBin');
end


%getting right lfp
for dat_cnt=8
    date=dates(dat_cnt,:);
    date
    all_data=cell(1,4);
    freq_list=cell(4,1);
    freq_list{1}='0.25';
    freq_list{2}='0.5';
    freq_list{3}='1';
    freq_list{4}='2';
    for ChN=1:60
        load(['F:\',date,'\processing\preprocessed\CH',int2str(ChN)]);
        for m=1:4
            for cond=1:3
                switch cond
                    case 1
                        cond_start=120;
                    case 2
                        cond_start=140;
                    case 3
                        cond_start=160;
                end
                freq=freq_list{m};
                switch freq
                    case '0.25'
                        per=2500;
                        start=1;
                    case '0.5'
                        per=1250;
                        start=2;
                    case '1'
                        per=625;
                        start=3;
                    case '2'
                        per=300;
                        start=4;
                end
                omit=20;
                load(strcat('D:\Adapt_experiment\stim_times_half3\hf3_',freq,'.mat'));
                stim=25*round(stim/omit);
                aver=mean(lfp(:,start+cond_start:4:cond_start+16),2);
                tmp=[];
                for i=3:2:length(stim)-1
                    tmp=[tmp aver(stim(i):stim(i)+per*2)];
                end
                tmp=mean(tmp,2);
                all_data{m}(cond,1:numel(tmp),ChN)=tmp;
            end
        end
    end
    descr=['for ',date,'; ND4(2) - 4 repetitions (without the last one). cell array by frequency; columns - time(variable); rows - condition; slices - channels (all 16)'];
    save(['F:\',date,'\processing\preprocessed\all_',date,'_cor'],'all_data','descr');
end

%getting right spikes
for dat_cnt=1:10
    date=dates(dat_cnt,:);
    date
    all_data=cell(1,4);
    freq_list=cell(4,1);
    freq_list{1}='0.25';
    freq_list{2}='0.5';
    freq_list{3}='1';
    freq_list{4}='2';
    shift=20;%used shift in ms
    for ChN=1:60
        load(['F:\',date,'\processing\preprocessed\CH',int2str(ChN),'_spikes']);
        for m=1:4
            for cond=1:3
                switch cond
                    case 1
                        cond_start=120;
                    case 2
                        cond_start=140;
                    case 3
                        cond_start=160;
                end
                freq=freq_list{m};
                switch freq
                    case '0.25'
                        per=100;
                        start=1;
                    case '0.5'
                        per=50;
                        start=2;
                    case '1'
                        per=25;
                        start=3;
                    case '2'
                        per=12;
                        start=4;
                end
                load(strcat('D:\Adapt_experiment\stim_times_half3\hf3_',freq,'.mat'));
                stim=floor(stim/shift);
                aver=mean(spikesBin(:,start+cond_start:4:cond_start+16),2);
                tmp=[];
                for i=3:2:length(stim)-1
                    tmp=[tmp aver(stim(i):stim(i)+per*2)];
                end
                tmp=mean(tmp,2);
                all_data{m}(cond,1:numel(tmp),ChN)=tmp;
            end
        end
    end
    descr=['for ',date,'; ND4(2) - 4 repetitions (without the last one). cell array by frequency; columns - time(variable); rows - condition; slices - channels (all 16)'];
    save(['F:\',date,'\processing\preprocessed\all_',date,'_spikes_cor'],'all_data','descr');
end