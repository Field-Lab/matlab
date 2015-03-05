clear
date='20121023';
unitspath=['/mnt/muench_data/user/alexandra/MEA_data/',date,'/easy_formatted_units/'];
units=dir([unitspath,'*FFFlicker.mat']);

path2save=['/mnt/muench_data/user/alexandra/MEA_data/',date,'/FFFlicker_LF_raw/'];

if ~exist(path2save,'file')
    mkdir(path2save);
end
load([unitspath,units(1).name])
filter_length=500; %filter length in ms
lseq=[];hseq=[];%m=[];
for i=1:length(spike_info.name_info)
    if ~isempty(cell2mat(regexp(spike_info.name_info(i),'LCseq')))
        lseq=[lseq i];
    else
        hseq=[hseq i];
    end
end
units_cnt=9;
units(units_cnt).name
load([unitspath,units(units_cnt).name])
trials=size(spike_info.flip_times,1);
% trials=144;

HighSpikeCount=zeros(1,length(hseq));
LowSpikeCount=zeros(1,length(lseq));
HighFilterRow=cell(length(hseq),1);
LowFilterRow=cell(length(lseq),1);
cntL=1;
cntH=1;

con=zeros(12,62000);
cnt=1;
for i=49:2:72
    spikes=spike_info.spike_times{i}'; % spikes in ms, rounded
%     spikes(spikes<filter_length|spikes>flips(end))=[]; 
    tmp=convolved(spikes,40,62000);
    con(cnt,:)=tmp(121:end-120);
    cnt=cnt+1;
end
 plot(con(1:3:end,:)')
imagesc(corr(con'))
[c,lags]=xcorr(con(1,:)',con(4,:)',500,'coeff');
[c,lags]=xcorr(con(3,:),500,'coeff');
plot(lags,c)

a=con(1,:);
thr=mean(a)+std(a);
b=a>thr;
c=diff(b);
spikes=find(c==1);

a=a-mean(a);
b=0;
for i=1:length(a)-500
    b=b+a(i+500)*tmp(i+499:-1:i);
end
b=b/(length(a)-500);
plot(b)


for trial_counter=1:trials
    flips=spike_info.flip_times{trial_counter}(:,1); % flips in ms, rounded
    pxls=(((spike_info.flip_times{trial_counter}(:,2))-30)/30); % pxs in ms, normalized -1:1
    spikes=spike_info.spike_times{trial_counter}'; % spikes in ms, rounded
    spikes(spikes<filter_length|spikes>flips(end))=[];    

    if isempty(spikes)
        if ~isempty(find(hseq==trial_counter, 1))
            HighSpikeCount(cntH)=0;
            cntH=cntH+1;            
        else
            LowSpikeCount(cntL)=0;
            cntL=cntL+1;
        end
    else
        
        tmp=zeros(flips(end)+1,2); % 1 column - color, 2 column - spikes
        tmp([1; flips(1:end-1)],1)=pxls;
        
        % fill stimulus down to ms with brightness values where spikes
        % happen
        for i=1:50
            subscr=(flips(1:end-1)+i);
            subscr(tmp(flips(1:end-1)+i,1)~=0)=[];
            tmp(subscr,1)=tmp(subscr-1,1);
        end
        
        tmp(spikes,2)=1;
        
        n=zeros(length(spikes),filter_length); % filter
        for i=1:filter_length
            n(:,i)=tmp(spikes-i+1,1);
            %     n(i)=sum(tmp(c-i+1,1))/all_spikes;
        end
        n=[spikes', n];
        % plot(b)
        tmp=n(n(:,1)<flips(3600)&n(:,1)>filter_length,2:end);
        if ~isempty(find(hseq==trial_counter, 1))
            HighFilterRow{cntH}=tmp;
            HighSpikeCount(cntH)=size(tmp,1);
            cntH=cntH+1;
        else
            LowFilterRow{cntL}=tmp;
            LowSpikeCount(cntL)=size(tmp,1);
            cntL=cntL+1;
        end
    end
end
save([path2save,units(units_cnt).name(1:end-27),'_FFFlicker_linear_filter_raw_double_spikes.mat'],'HighFilterRow','HighSpikeCount','LowFilterRow','LowSpikeCount')


HighSpikeCount=zeros(1,length(hseq));
LowSpikeCount=zeros(1,length(lseq));
HighFilterConv=zeros(700,length(hseq));
LowFilterConv=zeros(700,length(lseq));
cntL=1;
cntH=1;
sig=20;
for trial_counter=1:trials
    trial_counter
    flips=spike_info.flip_times{trial_counter}(:,1); % flips in ms, rounded
    pxls=(((spike_info.flip_times{trial_counter}(:,2))-30)/30); % pxs in ms, normalized -1:1
    spikes=spike_info.spike_times{trial_counter}'; % spikes in ms, rounded
    spikes(spikes>flips(end))=[];
    
    if isempty(spikes)
        a=zeros(1,63000);
    else
        a=convolved(spikes,sig,63000);
        a=a(sig*3+1:end-sig*3);
        a=a-mean(a);
    end
    tmp=zeros(flips(end)+1,1); % 1 column - color, 2 column - spikes
    tmp([1; flips(1:end-1)])=pxls;
    
    % fill stimulus down to ms with brightness values where spikes
    % happen
    for i=1:50
        subscr=(flips(1:end-1)+i);
        subscr(tmp(flips(1:end-1)+i,1)~=0)=[];
        tmp(subscr,1)=tmp(subscr-1,1);
    end
    tic
    b=0;
    for i=202:length(a)-700
        b=b+a(i+500)*tmp(i+699:-1:i);
    end
    b=b/(length(a)-901);
    toc
    if ~isempty(find(hseq==trial_counter, 1))
        HighFilterConv(:,cntH)=b;
        HighSpikeCount(cntH)=size(spikes,1);
        cntH=cntH+1;
    else
        LowFilterConv(:,cntL)=b;
        LowSpikeCount(cntL)=size(spikes,1);
        cntL=cntL+1;
    end
    
end
save([path2save,units(units_cnt).name(1:end-27),'_FFFlicker_linear_filter_conv_sig1.mat'],'HighFilterConv','HighSpikeCount','LowFilterConv','LowSpikeCount')


nd='87654321234567';

figure
set(gcf,'Position', [1 31 1680 946])
[rows,cols]=opt_subplots(length(nd));

minn=min(min(HighFilterConv(:,1:2:end)));
maxx=max(max(HighFilterConv(:,1:2:end)));
cnt=1;
for i=1:12:12*length(nd)
    k=subplot(rows,cols,cnt);
    hold on
    goto=min(i+11,size(HighFilterConv,2));
    plot(HighFilterConv(:,i:2:goto))
    hleg=legend(int2str(HighSpikeCount(i:2:goto)'));
    legend('boxoff')
    set(hleg,'Fontsize',8)
    title(['ND',nd(cnt)])
    axis([0 700 minn maxx])
    line([200 200],[minn maxx],'color','k')    
    line([0 700],[0 0],'color','k')
    cnt=cnt+1;
end





nd='87654321234567';

figure
set(gcf,'Position', [1 31 1680 946])
[rows,cols]=opt_subplots(12);

minn=min(min(HighFilterConv(:,1:2:end)));
maxx=max(max(HighFilterConv(:,1:2:end)));
cnt=1;
for i=61:72
    k=subplot(rows,cols,cnt);
    hold on
    a=HighFilterConv(:,i)/max(HighFilterConv(:,i));
    plot(a)
   a=LowFilterConv(:,i)/max(LowFilterConv(:,i));
    plot(a,'r')
    axis tight
    line([200 200],[minn maxx],'color','k')    
    line([0 700],[0 0],'color','k')
    cnt=cnt+1;
end


nd='87654321234567';

figure
set(gcf,'Position', [1 31 1680 946])
[rows,cols]=opt_subplots(length(nd));

minn=min(min(LowFilterConv(:,1:2:end)));
maxx=max(max(LowFilterConv(:,1:2:end)));
cnt=1;
for i=1:12:12*length(nd)
    k=subplot(rows,cols,cnt);
    hold on
    goto=min(i+11,size(LowFilterConv,2));
    plot(LowFilterConv(:,i:2:goto))
    hleg=legend(int2str(LowSpikeCount(i:2:goto)'));
    legend('boxoff')
    set(hleg,'Fontsize',8)
    title(['ND',nd(cnt)])
    axis([0 700 minn maxx])
    line([200 200],[minn maxx],'color','k')    
    line([0 700],[0 0],'color','k')
    cnt=cnt+1;
end
