cd('/mnt/muench_data/user/alexandra/scripts')

%% Figure M1
clear
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all')

figure
for i=1:8
    a(i)=mean(zc_HC(i,zc_HC(i,:)>0&onOff>0));
    b(i)=std(zc_HC(i,zc_HC(i,:)>0))/sqrt(sum(zc_HC(i,:)>0&onOff>0));
    c(i)=mean(zc_HC(i,zc_HC(i,:)>0&onOff<0));
    d(i)=std(zc_HC(i,zc_HC(i,:)>0))/sqrt(sum(zc_HC(i,:)>0&onOff<0));
end
errorbar(a,b,'r','linewidth',2)
hold on
errorbar(c,d,'b','linewidth',2)

a=peak_HC-peak_LC;
peak_HC(abs(a)>30)=0;
peak_LC(abs(a)>30)=0;
clear a

for i=1:8
    a(i)=mean(peak_HC(i,peak_HC(i,:)>0&onOff>0));
    b(i)=std(peak_HC(i,peak_HC(i,:)>0))/sqrt(sum(peak_HC(i,:)>0&onOff>0));
    c(i)=mean(peak_HC(i,peak_HC(i,:)>0&onOff<0));
    d(i)=std(peak_HC(i,peak_HC(i,:)>0))/sqrt(sum(peak_HC(i,:)>0&onOff<0));
end
errorbar(a,b,'m','linewidth',2)
hold on
errorbar(c,d,'c','linewidth',2)

legend({'zero crossing ON','zero crossing OFF','peak ON','peak OFF'})
axis([0 9 0 220])
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
title('Fig.M1. Latencies of peak and zero crossing for ON and OFF cells, mean +/-st error. High contrast GWN')


%% Fig.M2 temporal development of latency
clear
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all')
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/temporal_develop_LF.mat','HCon','HCoff','names','onOff','exp_codes')
[~, tmp]=min(HCoff);
tmp=reshape(tmp, size(tmp,2),size(tmp,3));
tmp(tmp<50)=nan;
tmp(tmp>250)=nan;
peak_off=nanmean(tmp');
peak_off_std=nanstd(tmp')./sqrt(sum(~isnan(tmp')));
for i=1:size(HCoff,2)
    for j=1:size(HCoff,3)
        if ~isnan(tmp(i,j))
            tmp1=HCoff(:,i,j);
            zcOFF(i,j)=find(tmp1(tmp(i,j):end)>0,1)+tmp(i,j)-1;
        else
            zcOFF(i,j)=nan;
        end
    end
end
zc_off=nanmean(zcOFF');
zc_off_std=nanstd(zcOFF')./sqrt(sum(~isnan(zcOFF')));


[~, tmp]=max(HCon);
tmp=reshape(tmp, size(tmp,2),size(tmp,3));
tmp(tmp<50)=nan;
tmp(tmp>250)=nan;
peak_on=nanmean(tmp');
peak_on_std=nanstd(tmp')./sqrt(sum(~isnan(tmp')));
for i=1:size(HCon,2)
    for j=1:size(HCon,3)
        if ~isnan(tmp(i,j))
            tmp1=HCon(:,i,j);
            zcON(i,j)=find(tmp1(tmp(i,j):end)<0,1)+tmp(i,j)-1;
        else
            zcON(i,j)=nan;
        end
    end
end
zc_on=nanmean(zcON');
zc_on_std=nanstd(zcON')./sqrt(sum(~isnan(zcON')));

c=5;coord=[];
for i=1:8
    coord=[coord c:c+8];
    c=c+20;
end




c=1;
for i=1:size(HCoff,3)
    HCoff(:,:,i)=HCoff(:,:,i)./repmat(sum(abs(HCoff(:,:,i))),500,1);
end
for i=1:8
    [~,peak_off(i,:)]=min(nanmean(HCoff(:,c:c+8,:),3));
    [~,peak_on(i,:)]=max(nanmean(HCon(:,c:c+8,:),3));
    for m=1:9
        zc_off(i,m)=find(nanmean(HCoff(mean(peak_off(i,:)):end,c+m-1,:),3)>0,1)+mean(peak_off(i,:))-1;
        zc_on(i,m)=find(nanmean(HCon(mean(peak_on(i,:)):end,c+m-1,:),3)<0,1)+mean(peak_on(i,:))-1;
    end
    c=c+9;
end


c=5;coord=[];
for i=1:8
    coord=[coord c:c+8];
    c=c+20;
end

errorbar(coord,peak_on,peak_on_std,'r.')
hold on
errorbar(coord,peak_off,peak_off_std,'b.')
errorbar(coord,zc_on,zc_on_std,'m.')
errorbar(coord,zc_off,zc_off_std,'c.')

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_M2','coord','peak_on','peak_off',...
    'peak_on_std','peak_off_std','zc_on','zc_off','zc_on_std','zc_off_std')

legend({'peak ON','peak OFF','zero crossing ON','zero crossing OFF'},'fontsize',12,'fontweight','bold')
axis([0 160 70 230])
set(gca,'xtick',9.5:20:160,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND','fontweight','bold','fontsize',12)
ylabel('Latency, ms','fontweight','bold','fontsize',12)

title('Fig.M2. Latencies of peak and zero crossing for ON and OFF cells in temporal development. High contrast GWN')


%% Fig.M5 start at nd4
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/temporal_develop_LF_nd4.mat','HCon','HCoff','names','onOff','exp_codes')
[~, tmp]=min(HCoff);
tmp=reshape(tmp, size(tmp,2),size(tmp,3));
tmp(tmp<50)=nan;
tmp(tmp>250)=nan;
peak_off=nanmean(tmp');
peak_off_std=nanstd(tmp')./sqrt(sum(~isnan(tmp')));


[~, tmp1]=max(HCon);
tmp1=reshape(tmp1, size(tmp1,2),size(tmp1,3));
tmp1(tmp1<50)=nan;
tmp1(tmp1>250)=nan;

peaks=nanmean([tmp'; tmp1']);
peaks_std=nanstd([tmp'; tmp1'])./sqrt(sum(~isnan([tmp'; tmp1'])));

c=3;coord=[];
for i=1:3
    coord=[coord c:c+3];
    c=c+7;
end


save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_M5','coord','peaks','peaks_std')

%% Fig.M5_alt start at ND4 (full version)

load('/mnt/muench_data/data/alexandra/MEA_data/analysis/temporal_develop_LF_nd4.mat','HCon','HCoff','names','onOff','exp_codes')
[~, tmp]=min(HCoff);
tmp=reshape(tmp, size(tmp,2),size(tmp,3));
tmp(tmp<50)=nan;
tmp(tmp>250)=nan;
peak_off=nanmean(tmp');
peak_off_std=nanstd(tmp')./sqrt(sum(~isnan(tmp')));
for i=1:size(HCoff,2)
    for j=1:size(HCoff,3)
        if ~isnan(tmp(i,j))
            tmp1=HCoff(:,i,j);
            zcOFF(i,j)=find(tmp1(tmp(i,j):end)>0,1)+tmp(i,j)-1;
        else
            zcOFF(i,j)=nan;
        end
    end
end
zc_off=nanmean(zcOFF');
zc_off_std=nanstd(zcOFF')./sqrt(sum(~isnan(zcOFF')));


[~, tmp]=max(HCon);
tmp=reshape(tmp, size(tmp,2),size(tmp,3));
tmp(tmp<50)=nan;
tmp(tmp>250)=nan;
peak_on=nanmean(tmp');
peak_on_std=nanstd(tmp')./sqrt(sum(~isnan(tmp')));
for i=1:size(HCon,2)
    for j=1:size(HCon,3)
        if ~isnan(tmp(i,j))
            tmp1=HCon(:,i,j);
            zcON(i,j)=find(tmp1(tmp(i,j):end)<0,1)+tmp(i,j)-1;
        else
            zcON(i,j)=nan;
        end
    end
end
zc_on=nanmean(zcON');
zc_on_std=nanstd(zcON')./sqrt(sum(~isnan(zcON')));

c=3;coord=[];
for i=1:3
    coord=[coord c:c+3];
    c=c+7;
end


save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_M5','coord','peak_on','peak_off',...
    'peak_on_std','peak_off_std','zc_on','zc_off','zc_on_std','zc_off_std')

%% Fig.M6 skip nd3
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/temporal_develop_LF_skipND3.mat','HCon','HCoff','names','onOff','exp_codes')
[~, tmp]=min(HCoff);
tmp=reshape(tmp, size(tmp,2),size(tmp,3));
tmp(tmp<50)=nan;
tmp(tmp>250)=nan;
peak_off=nanmean(tmp');
peak_off_std=nanstd(tmp')./sqrt(sum(~isnan(tmp')));
for i=1:size(HCoff,2)
    for j=1:size(HCoff,3)
        if ~isnan(tmp(i,j))
            tmp1=HCoff(:,i,j);
            zcOFF(i,j)=find(tmp1(tmp(i,j):end)>0,1)+tmp(i,j)-1;
        else
            zcOFF(i,j)=nan;
        end
    end
end
zc_off=nanmean(zcOFF');
zc_off_std=nanstd(zcOFF')./sqrt(sum(~isnan(zcOFF')));


[~, tmp]=max(HCon);
tmp=reshape(tmp, size(tmp,2),size(tmp,3));
tmp(tmp<50)=nan;
tmp(tmp>250)=nan;
peak_on=nanmean(tmp');
peak_on_std=nanstd(tmp')./sqrt(sum(~isnan(tmp')));
for i=1:size(HCon,2)
    for j=1:size(HCon,3)
        if ~isnan(tmp(i,j))
            tmp1=HCon(:,i,j);
            zcON(i,j)=find(tmp1(tmp(i,j):end)<0,1)+tmp(i,j)-1;
        else
            zcON(i,j)=nan;
        end
    end
end
zc_on=nanmean(zcON');
zc_on_std=nanstd(zcON')./sqrt(sum(~isnan(zcON')));

c=5;coord=[];
for i=1:5
    coord=[coord c:c+8];
    c=c+20;
end


save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_M6','coord','peak_on','peak_off',...
    'peak_on_std','peak_off_std','zc_on','zc_off','zc_on_std','zc_off_std')



%% Fig.M3 Gain evolution

clear a_on a_off
for i=1:8
    b=hc_gain(bothContr&onOff>0,i);  
    ons=ones(size(b))./b;
    ons=ons(ons>0&ons<35);
    a_on(i,1)=nanmean(ons);
    a_on(i,2)=nanstd(ons)/sqrt(length(ons));
    hold on
    
    b=hc_gain(bothContr&onOff<0,i);
    ons=ones(size(b))./b;
    ons=ons(ons>0&ons<35);
    a_off(i,1)=nanmean(ons);
    a_off(i,2)=nanstd(ons)/sqrt(length(ons));

end
figure
errorbar(a_on(:,1),a_on(:,2),'r')
hold on
errorbar(a_off(:,1),a_off(:,2))
legend({'ON','OFF'})
axis([0 9 5 25])
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND')
ylabel('gain, a.u.')
title('Fig.M3. Gain. High contrast GWN')

%% Fig.M4 Temporal development of gain
clear
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/temporal_develop_gain','nonLinear','names','onOff','exp_codes')

for i=1:72
    tmp=nonLinear(4,i,onOff>0);
    tmp=ones(size(tmp))./tmp;
    tmp=tmp(tmp>0&tmp<20);
    k(i,1)=nanmean(tmp);
    k(i,2)=nanstd(tmp)/sqrt(length(tmp));
    
    tmp=nonLinear(4,i,onOff<0);
    tmp=ones(size(tmp))./tmp;
    tmp=tmp(tmp>0&tmp<20);
    m(i,1)=nanmean(tmp);
    m(i,2)=nanstd(tmp)/sqrt(length(tmp));
end

figure
j=1;
for i=1:9:72
    errorbar(j:j+8,k(i:i+8,1),k(i:i+8,2),'r')
    hold on
    errorbar(j:j+8,m(i:i+8,1),m(i:i+8,2))
    j=j+12;
end
axis([0 98 0 20])
set(gca,'xtick',5:12:100,'xticklabel',{'8','7','6','5','4','3','2','1'})
legend('ON','OFF')
xlabel('ND')
ylabel('gain, a.u.')
title('Fig.M4. Temporal evolution of gain at high contrast GWN')

%% Fig.N

clear
date='20121023'
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord],'LinearFilter','SpikeCount','names')

figure
c=25;
cm=1;
data_acc=zeros(300,12,3);
for i=2:7
    subplot(3,6,i-1)
    a=reshape(LinearFilter(:,1,c+16:2:c+23,12),500,4);
    a=mean(a');
    a=a./max((a));    
    data_acc(:,cm,1)=a(1:300);
    plot(a,'b','linewidth',2)
    hold on
    a=reshape(LinearFilter(:,1,c+17:2:c+23,12),500,4);
    a=mean(a');
    a=a./max((a));
    data_acc(:,cm+1,1)=a(1:300);
    plot(a,'r','linewidth',2)
    legend('HC','LC')%,'Location','SouthEast')
    title(['ND', int2str(9-i)])
    line([0 500],[0 0],'color','k')
    axis([0 300 -1.5 1.5])

    subplot(3,6,i+5)
    a=reshape(LinearFilter(:,1,c+16:2:c+23,3),500,4);
    a=mean(a');
    a=a./max((a));
    data_acc(:,cm,2)=a(1:300);
    plot(a,'b','linewidth',2)
    hold on
    a=reshape(LinearFilter(:,1,c+17:2:c+23,3),500,4);
    a=mean(a');
    a=a./max((a));
    data_acc(:,cm+1,2)=a(1:300);
    plot(a,'r','linewidth',2)
    legend('HC','LC')%,'Location','SouthEast')
    title(['ND', int2str(9-i)])
    line([0 500],[0 0],'color','k')
    axis([0 300 -2.5 1.5])

    subplot(3,6,i+11)
    a=reshape(LinearFilter(:,1,c+16:2:c+23,13),500,4);
    a=mean(a');
    a=a./max((a));
    data_acc(:,cm,3)=a(1:300);
    plot(a,'b','linewidth',2)
    hold on
    a=reshape(LinearFilter(:,1,c+17:2:c+23,13),500,4);
    a=mean(a');
    a=a./max((a));
    data_acc(:,cm+1,3)=a(1:300);
    plot(a,'r','linewidth',2)
    legend('HC','LC')%,'Location','SouthEast')
    title(['ND', int2str(9-i)])
    line([0 500],[0 0],'color','k')
    c=c+24;
    axis([0 300 -1.5 1.5])
    cm=cm+2;
end

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_N.mat','data_acc')

%% Fig.N1
cd('/mnt/muench_data/user/alexandra/scripts')
clear
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summaryzc_late')
% load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all')
boundaries=[155 130 115 85 75 75 85 85];
figure
th=1.5;
cuts=20;
tmp=zc_HC-zc_LC;
for i=2:7
    subplot(2,6,i-1)
    cond=zc_LC(i,:)>0&zc_HC(i,:)>0&bothContr&zc_LC(i,:)~=boundaries(i)&zc_HC(i,:)~=boundaries(i);
    
    a=zc_HC(i,cond&onOff<0)-zc_LC(i,cond&onOff<0);
    a(abs(a)>cuts)=[];
    n(1,1)=sum(a<=-th)/length(a)*100;
    n(1,2)=sum(abs(a)<th)/length(a)*100;
    n(1,3)=sum(a>=th)/length(a)*100;
    
    a=zc_HC(i,cond&onOff>0)-zc_LC(i,cond&onOff>0);
    a(abs(a)>cuts)=[];
    n(2,1)=sum(a<=-th)/length(a)*100;
    n(2,2)=sum(abs(a)<th)/length(a)*100;
    n(2,3)=sum(a>=th)/length(a)*100;
    
    bar(n')
    axis([0 4 0 100])
    title(['ND',int2str(9-i),'  zero crossing'])
    
    subplot(2,6,i+5)
    a=peak_HC(i,cond&onOff<0)-peak_LC(i,cond&onOff<0);
    a(abs(a)>cuts)=[];
    n(1,1)=sum(a<=-th)/length(a)*100;
    n(1,2)=sum(abs(a)<th)/length(a)*100;
    n(1,3)=sum(a>=th)/length(a)*100;
    
    a=peak_HC(i,cond&onOff>0)-peak_LC(i,cond&onOff>0);
    a(abs(a)>cuts)=[];
    n(2,1)=sum(a<=-th)/length(a)*100;
    n(2,2)=sum(abs(a)<th)/length(a)*100;
    n(2,3)=sum(a>=th)/length(a)*100;
    bar(n')
    title(['ND',int2str(9-i),'  peak'])
    axis([0 4 0 100])
end

%% Fig.N1 alternative - distribution
cd('/mnt/muench_data/user/alexandra/scripts')
clear
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summaryzc_late')
% load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all')
boundaries=[155 130 115 85 75 75 85 85];
figure
cuts=30;
for i=1:8
    subplot(4,4,i)
    cond=zc_LC(i,:)>0&zc_HC(i,:)>0&bothContr&zc_LC(i,:)~=boundaries(i)&zc_HC(i,:)~=boundaries(i);
    
    a=zc_HC(i,cond&onOff<0)-zc_LC(i,cond&onOff<0);
    a(abs(a)>cuts)=[];
    [m,k]=hist([a -40 40],13);
    plot(k(2:end-1),m(2:end-1)/length(a)*100,'b')
    hold on
    a=zc_HC(i,cond&onOff>0)-zc_LC(i,cond&onOff>0);
    a(abs(a)>cuts)=[];
    [m,k]=hist([a -40 40],13);
    plot(k(2:end-1),m(2:end-1)/length(a)*100,'r')
    title(['ND',int2str(9-i),' zc'])
    legend('OFF','ON')
    axis([-30 30 0 60])
    line([0,0],[0,60],'color','k')
    
    subplot(4,4,i+8)
    hold on
    a=peak_HC(i,cond&onOff<0)-peak_LC(i,cond&onOff<0);
    a(abs(a)>cuts)=[];
    [m,k]=hist([a -40 40],13);
    plot(k(2:end-1),m(2:end-1)/length(a)*100)
    a=peak_HC(i,cond&onOff>0)-peak_LC(i,cond&onOff>0);
    a(abs(a)>cuts)=[];
    [m,k]=hist([a -40 40],13);
    plot(k(2:end-1),m(2:end-1)/length(a)*100,'r')
    title(['ND',int2str(9-i),' peak'])
    legend('OFF','ON')
    ylabel('%')
    xlabel('HC-LC')
    axis([-30 30 0 60])
    line([0,0],[0,60],'color','k')
end

%% Fig.N1 alt alt
clear
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all')
for i=1:8
    for j=find(bothContr)
        [ind,lags]=xcorr(HC(50:300,i,j),LC(50:300,i,j),'coeff');
        [mc,ic]=max(ind);
        delay(j,i)=lags(ic);
        value(j,i)=mc;
    end
end

figure
a_on=cell(8,1);
a_off=cell(8,1);
data_acc_on=zeros(75,2,8);
data_acc_off=zeros(75,2,8);
for i=1:8
    subplot(2,4,i)
    a=delay(onOff>0,i);
    b=value(onOff>0,i);
    c=a<40&a>-40&b>0.6;
    [k,kind]=hist([a(c)' -100 100],75);
    plot(kind,k/sum(c)*100,'-*r')
    data_acc_on(:,1,i)=k/sum(c)*100;
    data_acc_on(:,2,i)=kind;
    
    
    axis([-40 20 0 40])    
    a_on{i}=a(c);
    hold on
    a=delay(onOff<0,i);
    b=value(onOff<0,i);
    c=a<40&a>-40&b>0.6;
    [k1,kind]=hist([a(c)' -100 100],75);
    plot(kind,k1/sum(c)*100,'-*')
    data_acc_off(:,1,i)=k1/sum(c)*100;
    data_acc_off(:,2,i)=kind;
    
    a_off{i}=a(c);
    ranksum(a_on{i},a_off{i})
    line([0,0],[0 40],'color','k')
    title(['ND',int2str(9-i)])
    xlabel('Delay of max crosscorr')
    ylabel('%')

end

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_N_alt_alt.mat','data_acc_on','data_acc_off')


clear p_on p_off
for i=1:7
    for j=i+1:8
        p_on(i,j-1) = ranksum(a_on{i},a_on{j});
        p_off(i,j-1) = ranksum(a_off{i},a_off{j});
    end
end

figure
subplot(1,2,1)
a=p_on;
a(a>=0.05)=2;
a(a>0&a<0.05)=1;
imagesc(a)
set(gca,'xtick',1:7,'xticklabel',{'7','6','5','4','3','2','1'})
set(gca,'ytick',1:7,'yticklabel',{'8','7','6','5','4','3','2'})

subplot(1,2,2)
a=p_off;
a(a>=0.05)=2;
a(a>0&a<0.05)=1;
imagesc(a)
set(gca,'xtick',1:7,'xticklabel',{'7','6','5','4','3','2','1'})
set(gca,'ytick',1:7,'yticklabel',{'8','7','6','5','4','3','2'})

%% Fig.N2
cd('/mnt/muench_data/user/alexandra/scripts')
clear
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summaryzc_late')
% load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all')
boundaries=[155 130 115 85 75 75 85 85];
figure
th=1.5;
cuts=20;
tmp=zc_HC-zc_LC;
for i=2:7
    cond=zc_LC(i,:)>0&zc_HC(i,:)>0&bothContr&zc_LC(i,:)~=boundaries(i)&zc_HC(i,:)~=boundaries(i);
    
    a=zc_HC(i,cond&onOff<0)-zc_LC(i,cond&onOff<0);
    a(abs(a)>cuts)=[];
    a=a(a<=-th);
    n(i-1,1)=mean(a);
    n(i-1,2)=std(a)/sqrt(length(a));

    a=zc_HC(i,cond&onOff>0)-zc_LC(i,cond&onOff>0);
    a(abs(a)>cuts)=[];
    a=a(a<=-th);
    n1(i-1,1)=mean(a);
    n1(i-1,2)=std(a)/sqrt(length(a));
end
errorbar(n(:,1),n(:,2))
hold on
errorbar(n1(:,1),n1(:,2),'r')
axis([0 7 -16 16])

for i=2:7
    cond=zc_LC(i,:)>0&zc_HC(i,:)>0&bothContr&zc_LC(i,:)~=boundaries(i)&zc_HC(i,:)~=boundaries(i);
    
    a=zc_HC(i,cond&onOff<0)-zc_LC(i,cond&onOff<0);
    a(abs(a)>cuts)=[];
    a=a(a>=th);
    n(i-1,1)=mean(a);
    n(i-1,2)=std(a)/sqrt(length(a));
length(a)

    a=zc_HC(i,cond&onOff>0)-zc_LC(i,cond&onOff>0);
    a(abs(a)>cuts)=[];
    a=a(a>=th);
    n1(i-1,1)=mean(a);
    n1(i-1,2)=std(a)/sqrt(length(a));
    
end
errorbar(n(:,1),n(:,2),'c')
errorbar(n1(:,1),n1(:,2),'m')
legend('reg OFF','reg ON','wrong OFF','wrong ON')
axis([0 9 -16 16])
set(gca,'xtick',1:6,'xticklabel',{'7','6','5','4','3','2'})
xlabel('ND')
line([0 9],[0 0],'color','k')
title('Fig N2. Adaptation of regular and wrong ON anf OFF cells by ND: zero crossing')


%% Fig.N3
cd('/mnt/muench_data/user/alexandra/scripts')
clear
% load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summaryzc_late')
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all')
boundaries=[155 130 115 85 75 75 85 85];
figure
th=1.5;
cuts=20;
tmp=zc_HC-zc_LC;
for i=2:7
    cond=zc_LC(i,:)>0&zc_HC(i,:)>0&bothContr&zc_LC(i,:)~=boundaries(i)&zc_HC(i,:)~=boundaries(i);
    subplot(2,6,i-1)
    hold on
    a=zc_HC(i,cond&onOff<0)-zc_LC(i,cond&onOff<0);
    a(abs(a)>cuts)=[];
    a=a(a<=-th);
    n(i-1,1)=mean(a);
    n(i-1,2)=std(a)/sqrt(length(a));
    plot(zc_HC(i,cond&onOff<0),'.');    
    plot(zc_LC(i,cond&onOff<0),'.r');
    line( [0 100 ],[mean(zc_HC(i,cond&onOff<0)) mean(zc_HC(i,cond&onOff<0))],'color','b')
    line( [0 100 ],[mean(zc_LC(i,cond&onOff<0)) mean(zc_LC(i,cond&onOff<0))],'color','r')
    axis([0 100 0 250])
    
    a=zc_HC(i,cond&onOff>0)-zc_LC(i,cond&onOff>0);
    a(abs(a)>cuts)=[];
    a=a(a<=-th);
    n1(i-1,1)=mean(a);
    n1(i-1,2)=std(a)/sqrt(length(a));
    k1(i-1)=mean(zc_HC(i,cond&onOff<0));
    subplot(2,6,i+5)
    hold on
    plot(zc_HC(i,cond&onOff>0),'.');
    plot(zc_LC(i,cond&onOff>0),'.r');
    line( [0 100 ],[mean(zc_HC(i,cond&onOff>0)) mean(zc_HC(i,cond&onOff>0))],'color','b')
    line( [0 100 ],[mean(zc_LC(i,cond&onOff>0)) mean(zc_LC(i,cond&onOff>0))],'color','r')
    axis([0 100 0 250])
end

figure
th=1.5;
cuts=20;
tmp=zc_HC-zc_LC;
for i=2:7
    cond=zc_LC(i,:)>0&zc_HC(i,:)>0&bothContr&zc_LC(i,:)~=boundaries(i)&zc_HC(i,:)~=boundaries(i);
    subplot(2,6,i-1)
    hold on
    a=zc_HC(i,cond&onOff<0);
    b=a;    
    
    cond=zc_LC(i+1,:)>0&zc_HC(i+1,:)>0&bothContr&zc_LC(i+1,:)~=boundaries(i+1)&zc_HC(i+1,:)~=boundaries(i+1);
    a=zc_LC(i+1,cond&onOff<0);
    plot(a,'.r')
    plot(b,'.b')
    
    
    subplot(2,6,i+5)
       cond=zc_LC(i,:)>0&zc_HC(i,:)>0&bothContr&zc_LC(i,:)~=boundaries(i)&zc_HC(i,:)~=boundaries(i);

    hold on
    a=zc_HC(i,cond&onOff>0);
    b=a;    
    
    cond=zc_LC(i+1,:)>0&zc_HC(i+1,:)>0&bothContr&zc_LC(i+1,:)~=boundaries(i+1)&zc_HC(i+1,:)~=boundaries(i+1);
    a=zc_LC(i+1,cond&onOff>0);
    plot(a,'.r')
    plot(b,'.b')
end

%% Fig N4 Gain changes
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_gain','hc_gain','lc_gain')
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all','onOff','bothContr')


clear a b ons offs
for i=1:8
    
    a=hc_gain(bothContr&onOff>0,i);
    b=lc_gain(bothContr&onOff>0,i);
    m=a./b;
    m=m(m<20&m>0);
    ons(i,1)=nanmean(m);
    ons(i,2)=nanstd(m)/sqrt(length(m));
    a_on{i}=m;
    
    a=hc_gain(bothContr&onOff<0,i);
    b=lc_gain(bothContr&onOff<0,i);
    m=a./b;
    m=m(m<20&m>0);
    offs(i,1)=nanmean(m);
    offs(i,2)=nanstd(m)/sqrt(length(m));
    a_off{i}=m;
end
figure
errorbar(offs(:,1),offs(:,2))
hold on
errorbar(ons(:,1),ons(:,2),'r')

clear p
for i=1:8
    p(i)=ranksum(a_on{i},a_off{i});
end
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_N4','ons','offs','p')


legend('ON','OFF','location','best')
axis([0 9 0 10])
line([0 9],[5,5],'color','k')
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND')
ylabel('gain ratio')
title('Fig N4. Gain increase at low contrast compared to high')


clear p_on p_off
for i=1:7
    for j=i+1:8
        p_on(i,j-1) = ranksum(a_on{i},a_on{j});
        p_off(i,j-1) = ranksum(a_off{i},a_off{j});
    end
end

%% Fig N6 Gain ratio distribution: evaluate fig N4 first (previous cell)

figure
data_acc=zeros(10,4,8);
for i=1:8
    subplot(2,4,i)
    tmp=a_on{i};
    [k,kind]=hist([tmp' 0 20],10);
    plot(kind,k/length(tmp)*100,'-*r')
    axis([0 20 0 50]) 
    data_acc(:,1,i)=kind;
    data_acc(:,2,i)=k/length(tmp)*100;

    hold on
    tmp=a_off{i};
    [k,kind]=hist([tmp' 0 20],10);
    plot(kind,k/length(tmp)*100,'-*')
    data_acc(:,3,i)=kind;
    data_acc(:,4,i)=k/length(tmp)*100;
    
    line([5,5],[0 50],'color','k')
    title(['ND',int2str(9-i)])
    xlabel('Gain ratio')
    ylabel('%')
end

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_N6','data_acc')




%% Fig N5 Delay vs Firing rate change

clear
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all')
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_gain','hc_gain','lc_gain')
for i=1:8
    for j=find(bothContr)
        [ind,lags]=xcorr(HC(50:300,i,j),LC(50:300,i,j),'coeff');
        [mc,ic]=max(ind);
        delay(j,i)=lags(ic);
        value(j,i)=mc;
        gainContr(j,i)=hc_gain(j,i)./lc_gain(j,i);
        fr(j,i)=frMean_HC(j,i)-frMean_LC(j,i);
    end
end

data=cell(2,8);
figure
for i=1:8
    subplot(2,4,i)
    a=delay(onOff>0,i);
    b=value(onOff>0,i);
    c=fr(onOff>0,i);
    d=a<30&a>-30&b>0.6;
    a=a(d);
%     a=a/max(abs(a));
    c=c(d);
%     c=c/max(abs(c));

    data{1,i}=[a c];
    plot(a,c,'or')
    hold on
    
    a=delay(onOff<0,i);
    b=value(onOff<0,i);
    c=fr(onOff<0,i);
    d=a<30&a>-30&b>0.6;
    a=a(d);
%     a=a/max(abs(a));
    c=c(d);
%     c=c/max(abs(c));
    
    data{2,i}=[a c];
    plot(a,c,'o')

%     line([-1 1],[0 0],'color','k')
%     line([0,0],[-1 1],'color','k')
%     axis([-1 1 -1 1])

    title(['ND',int2str(9-i)])
    xlabel('Delay')
    ylabel('firing change')

end

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_N5','data')


%% Fig N9 Delay in non-typical, typical, and off cells

clear
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all')
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_gain','hc_gain','lc_gain')
for i=1:8
    for j=find(bothContr)
        [ind,lags]=xcorr(HC(50:300,i,j),LC(50:300,i,j),'coeff');
        [mc,ic]=max(ind);
        delay(j,i)=lags(ic);
        value(j,i)=mc;
        gainContr(j,i)=hc_gain(j,i)./lc_gain(j,i);
        fr(j,i)=frMean_HC(j,i)-frMean_LC(j,i);
    end
end

data=cell(2,8);

for i=1:8
    a=delay(onOff>0,i);
    b=value(onOff>0,i);
    c=fr(onOff>0,i);
    d=a<30&a>-30&b>0.6;
    a=a(d);
    c=c(d);

    data{1,i}=[a c];

    
    a=delay(onOff<0,i);
    b=value(onOff<0,i);
    c=fr(onOff<0,i);
    d=a<30&a>-30&b>0.6;
    a=a(d);
    c=c(d);
    
    data{2,i}=[a c];

end

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_N9','data')


%% Fig.N10
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_gain','hc_gain','lc_gain')
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all','onOff','bothContr')


for i=1:8
    for j=find(bothContr)
        [ind,lags]=xcorr(HC(50:300,i,j),LC(50:300,i,j),'coeff');
        [mc,ic]=max(ind);
        delay(j,i)=lags(ic);
        value(j,i)=mc;
        gainContr(j,i)=hc_gain(j,i)./lc_gain(j,i);
        fr(j,i)=frMean_HC(j,i)-frMean_LC(j,i);
    end
end


clear a b ons offs
for i=1:8
    a=delay(bothContr&onOff>0,i);
    b=value(bothContr&onOff>0,i);
    c=fr(onOff>0,i);
    d=a<30&a>-30&b>0.6;
    c=c(d);
    
    a=hc_gain(bothContr&onOff>0,i);
    b=lc_gain(bothContr&onOff>0,i);
    m=a./b;
    m=m(d);
    res_t{i}=m(c>0);
    res_nt{i}=m(c<0);    
end

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_N10','res_t','res_nt')


%% Fig.Z11 Z12 Z13
for i=1:8
    
    tmp=bf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    a=std(tmp(50:450,:));
    for j=1:200
        [k(i,j),m(i,j)]=max(tmp(500:1000,j));
        b=find(tmp(500+m(i,j):2500,j)<a(j),1);
        if ~isempty(b)
            returns(j,i)=b+m(i,j);
        else
            returns(j,i)=2000;
        end
    end

end
ampl=k';
lat=m';

for i=1:8
    for j=1:8
        [~,p(i,j)]=ttest(ampl(:,i),ampl(:,j));
    end
end
clear a b
for i=1:8
    a(i)=mean(ampl(:,i));
    b(i)=std(ampl(:,i))/sqrt(200);
    
end
errorbar(a,b)

clear a b
for i=1:8
    a(i)=mean(lat(:,i));
    b(i)=std(lat(:,i))/sqrt(200);    
end
errorbar(a,b)
for i=1:8
    for j=1:8
        [~,p(i,j)]=ttest(lat(:,i),lat(:,j));
    end
end
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_Z11','lat','ampl','returns');


clear data
for i=1:8
    subplot(2,4,i)
    [k, kind]=hist([returns(:,i); 0; 2000],20);
    k(1)=k(1)-1;
    k(end)=k(end)-1;
    plot(kind,k/2,'-*');
    data(:,1,i)=kind;
    data(:,2,i)=k/2;
    axis([0 2000 0 50])
end

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_Z13','data');

% always sustained
for i=3:5
    kk(:,i-2)=returns(:,i)==2000;
end
sum(sum(kk(:,:)')==3)/2

% always transient
for i=3:5
    kk(:,i-2)=returns(:,i)<2000;
end
sum(sum(kk(:,:)')==3)/2

%% Fig.R13 
% fig Z6_alt
clear p h n sust_fav
figure
data=cell(8,1);
for i=1:8
    
    tmp=wf{i}(:,onOff>0);
    m=sum(tmp(500:1000,:))./sum(tmp(500:2500,:));    
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    [a lat]=max(tmp(500:1000,:));
    
    tmp=bf{i}(:,onOff>0);
    mOn=sum(tmp(2500:3000,:))./sum(tmp(2500:4500,:));
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    [aOn latOn]=max(tmp(2500:3000,:));
      
    %amplitude
    subplot(2,4,i)
    m=aOn./a;
    k=latOn-lat;
    b=isinf(m);
    data{i}=[k(~b);m(~b)];
    plot(k(~b),m(~b),'o')
    axis([-50 50 -0 2])
    xlabel('latency delta,ms. White-Black')
    ylabel('amplitude ratio. White/Black')
    title(['ND', int2str(9-i),])
    line([0,0],[0,2],'color','k')
    line([-50 50],[1,1],'color','k')

    p(i)=ranksum(a,aOn);
    h(i)=ranksum(lat,latOn);
    n(i)=ranksum(m,mOn);
 
end
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_R13','data')





%% Fig. Z
clear
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all')


% Fig Z
figure
nds='87654321'
data=zeros(4500,2,8);
for i=1:8
    subplot(2,4,i)
    tmp=wf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    hold off
    plot(mean(tmp,2),'color','k','linewidth',2)
    data(:,1,i)=mean(tmp,2); 
    
    hold on
    tmp=bf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    plot(4701:9200,mean(tmp,2),'color','k','linewidth',2)
    data(:,2,i)=mean(tmp,2);
 
    
    axis([1 9200 -5 65])
    line([1 4500],[0,0],'color','k')
    line([4700 9200],[0,0],'color','k')
    
    line([1 500],[55,55],'color','k')
    line([500 2500],[60,60],'color','k')
    line([2500 4500],[55,55],'color','k')
    line([500 500],[55,60],'color','k')
    line([2500 2500],[55,60],'color','k')
    
    line([1 500]+4700,[55,55],'color','k')
    line([500 2500]+4700,[50,50],'color','k')
    line([2500 4500]+4700,[55,55],'color','k')
    line([500 500]+4700,[50,55],'color','k')
    line([2500 2500]+4700,[50,55],'color','k')

    title(['ND',nds(i)])
end

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_Z','data')




nds='87654321';
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/OFFcellsClust.txt')
manual=zeros(520,1)-1;
rebSt=0;onSt=0;
on=[0 0 0];reb=[0 0 0];
rebStIs=0;onStIs=0;
feat=zeros(6,520);
for i=1:200
    b=OFFcellsClust(i,1);
    a=int2str(OFFcellsClust(i,2));
    while length(a)<6
        a=['0' a];        
    end
    if length(unique(a(4:6)))==1
        rebSt=rebSt+1;
        if unique(a(4:6))=='1'
            rebStIs=rebStIs+1;
        end
    end

    
    if length(unique(a(1:3)))==1
        onSt=onSt+1;
        if unique(a(1:3))=='1'
            onStIs=onStIs+1;
        end
    end
    
    if a(1)=='1'
        on(1)=on(1)+1;
        feat(1,b)=1;
    end
    if a(2)=='1'
        on(2)=on(2)+1;
        feat(2,b)=1;
    end
    if a(3)=='1'
        on(3)=on(3)+1;
        feat(3,b)=1;
    end
    if a(4)=='1'
        reb(1)=reb(1)+1;
        feat(4,b)=1;
    end
    if a(5)=='1'
        reb(2)=reb(2)+1;
        feat(5,b)=1;
    end
    if a(6)=='1'
        reb(3)=reb(3)+1;
        feat(6,b)=1;
    end
    manual(OFFcellsClust(i,1))=bin2dec(a);  
    
end




%Fig.Z1
figure
for i=1:8
    
    tmp=bf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    [a lat]=max(tmp(500:1000,:));    
    [aOn latOn]=max(tmp(2500:3000,:));
    
    subplot(2,4,i)
    plot(lat,a,'o')
    hold on
    plot(latOn, aOn,'or')

    axis([0 500 0 200])
    xlabel('latency,ms')
    ylabel('amplitude, Hz')
    title(['ND', int2str(9-i),])
    if i==1
        legend('OFF response','ON response')
    end
end



% Fig Z1_alt

figure
i=3
cc=1;
for cnt=[94 387 448]
    subplot(1,3,cc)
    plot(bf{i}(:,cnt),'linewidth',2)
    line([1 500],[140,140],'color','k')
    line([500 2500],[130,130],'color','k')
    line([2500 4500],[140,140],'color','k')
    line([500 500],[130,140],'color','k')
    line([2500 2500],[140,130],'color','k')
    axis([0 4500 0 150])
    set(gca,'xtick',1000:1000:4000,'xticklabel',{'1','2','3','4'})
    xlabel('time,s')
    ylabel('Firing rate, Hz')
    cc=cc+1;
end
data=bf{i}(:,[94 387 448]);
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_Z1_alt','data')

% Z2 Mean and std of the response to black flash, OFF cells
figure
for i=3:5
    subplot(1,3,i-2)
    tmp=bf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    b=mean(tmp');
    a=std(tmp');
    plot(b,'k','Linewidth',2)
    hold on
    plot(b-a,'k','Linewidth',1)
    plot(b+a,'k','Linewidth',1)
    title(['ND', nds(i)])
    axis([0 4500 -20 100])
    line([1 4500],[0,0],'color','k')
   
    line([1 500],[90,90],'color','k')
    line([500 2500],[85,85],'color','k')
    line([2500 4500],[90,90],'color','k')
    line([500 500],[85,90],'color','k')
    line([2500 2500],[85,90],'color','k')
    
end

% Z2_alt Std of the response to black flash, OFF cells
data=zeros(4500,3);
figure
for i=3:5
    subplot(1,3,i-2)
    tmp=bf{i}(:,onOff<0);
    a=std(tmp');
    data(:,i-2)=a;
    plot(a,'k','Linewidth',1)
    title(['ND', nds(i)])
    axis([0 4500 0 45])   
    line([1 500],[90,90]-50,'color','k')
    line([500 2500],[85,85]-50,'color','k')
    line([2500 4500],[90,90]-50,'color','k')
    line([500 500],[85,90]-50,'color','k')
    line([2500 2500],[85,90]-50,'color','k')
    set(gca,'xtick',1000:1000:4000,'xticklabel',{'1','2','3','4'})
    xlabel('time,s')
    ylabel('std')
end

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_Z2_alt','data')



% fig Z3
nds='87654321';
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/OFFcellsClust.txt')
manual=zeros(520,1)-1;
rebSt=0;onSt=0;
on=[0 0 0];reb=[0 0 0];
rebStIs=0;onStIs=0;
feat=zeros(6,520);
for i=1:200
    b=OFFcellsClust(i,1);
    a=int2str(OFFcellsClust(i,2));
    while length(a)<6
        a=['0' a];        
    end
    if length(unique(a(4:6)))==1
        rebSt=rebSt+1;
        if unique(a(4:6))=='1'
            rebStIs=rebStIs+1;
        end
    end

    
    if length(unique(a(1:3)))==1
        onSt=onSt+1;
        if unique(a(1:3))=='1'
            onStIs=onStIs+1;
        end
    end
    
    if a(1)=='1'
        on(1)=on(1)+1;
        feat(1,b)=1;
    end
    if a(2)=='1'
        on(2)=on(2)+1;
        feat(2,b)=1;
    end
    if a(3)=='1'
        on(3)=on(3)+1;
        feat(3,b)=1;
    end
    if a(4)=='1'
        reb(1)=reb(1)+1;
        feat(4,b)=1;
    end
    if a(5)=='1'
        reb(2)=reb(2)+1;
        feat(5,b)=1;
    end
    if a(6)=='1'
        reb(3)=reb(3)+1;
        feat(6,b)=1;
    end
    manual(OFFcellsClust(i,1))=bin2dec(a);  
    
end


a=manual;
a_list=[];
code_list=[];
for code=0:63
    if sum(a==code)>12
        code_list=[code_list code];
        a_list=[a_list; dec2bin(code,6)];
    end
end
figure
clear means
cc=1;tt=1;
data=zeros(4500,3,6);
for code=code_list([1 5 2 6 3 4])
    cc=tt;
    for i=3:5
        exx=int2str(length(unique(exp_codes(a==code))));
        subplot(3,6,cc)
        k=bf{i}';
        k=k-repmat(mean(k(:,50:450),2),1,4500);
        means(1+4500*(i-2):4500+4500*(i-2),tt)=mean(k(a==code,:));
        
        plot(mean(k(a==code,:)),'color','k','linewidth',2)
        
        data(:,i-2,tt)=mean(k(a==code,:));
        
        axis([0 4500 -10 75])
        hold on
        title(['ND',nds(i),' n=',int2str(sum(a==code)), ', ',exx,' exp'])
        
        line([1 500],[65,65],'color','k')
        line([500 2500],[60,60],'color','k')
        line([2500 4500],[65,65],'color','k')
        line([500 500],[60,65],'color','k')
        line([2500 2500],[65,60],'color','k')        
        
        set(gca,'xtick',1000:1000:4000,'xticklabel',{'1','2','3','4'})
        xlabel('time,s')
        ylabel('Firing Rate, Hz')
        
        cc=cc+6;
    end
    tt=tt+1;
end

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_Z3','data')


% Fig Z4
data=[on/2 reb/2];
figure
bar(on/2)
hold on
bar(6:8,reb/2,'facecolor',[0.5 0.1 0])
set(gca,'xtick',[1:3 6:8],'xticklabel',{'ND6','ND5','ND4','ND6','ND5','ND4'})
text(1,-5,'Early ON response')
text(6,-5,'Delayed ON response')
ylabel('% of OFF cells recorded')
title('Fig.Z4. Fraction of cells with early and delayed ON response')

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_Z4','data')

%Fig Z5
figure
a=[rebStIs rebSt-rebStIs]/2
b=[onStIs onSt-onStIs]/2
data=[b; a];
bar([b; a],'BarLayout','stacked')
set(gca,'xtick',1:2,'xticklabel',{'Early ON response','Delayed ON response'})
ylabel('% of OFF cells recorded') 
axis([0 3 0 70])
legend('always present','never present')
title('Fig.Z5. % of OFF cells with stably present or absent ON response')
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_Z5','data')

% fig Z6
clear p h n sust_fav
figure
for i=2:7
    
    tmp=bf{i}(:,onOff<0);
    m=sum(tmp(500:1000,:))./sum(tmp(550:2500,:));    
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    [a lat]=max(tmp(500:1000,:));
   
      
    %amplitude
    subplot(3,6,i-1)
    [k,kind]=hist([a -10 200],21);      
    plot(kind(2:end-1),k(2:end-1)/length(a)*100,'-*')
    hold on
    
    tmp=wf{i}(:,onOff<0);
    m1=sum(tmp(2500:3000,:))./sum(tmp(2550:4500,:));
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    [aOn latOn]=max(tmp(2500:3000,:));
    
    [k,kind]=hist([aOn -110 200],31);
    plot(kind(12:end-1),k(12:end-1)/length(aOn)*100,'-*r')
    
    p(i)=ranksum(a,aOn);
    h(i)=ranksum(lat,latOn);
    n(i)=ranksum(m,m1);
    
    axis([0 160 0 25])
    xlabel('amplitude, Hz')
    ylabel('%')
    title(['ND', int2str(9-i),])
    
    if i==7
        legend('black','white')
    end
    
    % latency
    subplot(3,6,i+5)
    [k,kind]=hist([lat 0 500],50); 
    k([1,end])=k([1,end])-1;
    k(k<0)=0;
    plot(kind,k/length(a)*100,'-*')
    hold on    
    [k,kind]=hist([latOn 0 500],50);  
    k([1,end])=k([1,end])-1;
    k(k<0)=0;
    plot(kind,k/length(a)*100,'-*r')    
    
    axis([70 300 0 30]) 
    xlabel('latency, ms')
    ylabel('%')
    title(['ND', int2str(9-i),])
    
    % transiency
    subplot(3,6,i+11)
    [k,kind]=hist([m 0 1],20); 
    k([1,end])=k([1,end])-1;
    k(k<0)=0;
    sust_fav(i,1)=sum(m<0.5)/length(m)*100;
    plot(kind,k/length(a)*100,'-*')
    hold on    
    [k,kind]=hist([m1 0 1],20); 
    k([1,end])=k([1,end])-1;
    k(k<0)=0;
    sust_fav(i,2)=sum(m1<0.5)/length(m1)*100;
    plot(kind,k/length(a)*100,'-*r')    
    
    axis([0 1 0 30]) 
    xlabel('transiency index')
    ylabel('%')
    title(['ND', int2str(9-i),])
end


% fig Z6_alt
clear p h n sust_fav
figure
data=cell(8,1);
for i=1:8
    
    tmp=bf{i}(:,onOff<0);
    m=sum(tmp(500:1000,:))./sum(tmp(500:2500,:));    
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    [a lat]=max(tmp(500:1000,:));
    
    tmp=wf{i}(:,onOff<0);
    mOn=sum(tmp(2500:3000,:))./sum(tmp(2500:4500,:));
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    [aOn latOn]=max(tmp(2500:3000,:));
      
    %amplitude
    subplot(2,4,i)
    m=aOn./a;
    k=latOn-lat;
    b=isinf(m);
    data{i}=[k(~b);m(~b)];
    plot(k(~b),m(~b),'o')
    axis([-50 50 -0 2])
    xlabel('latency delta,ms. White-Black')
    ylabel('amplitude ratio. White/Black')
    title(['ND', int2str(9-i),])
    line([0,0],[0,2],'color','k')
    line([-50 50],[1,1],'color','k')

    p(i)=ranksum(a,aOn);
    h(i)=ranksum(lat,latOn);
    n(i)=ranksum(m,mOn);
 
end
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_Z6_alt','data')


% fig Z7
figure
for i=1:8
    
    tmp=bf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    a=std(tmp(50:450,:));
    for j=1:200
        [k(i,j),m]=max(tmp(500:1000,j));
        b=find(tmp(500+m:2500,j)<a(j),1);
        if ~isempty(b)
            returns(j,i)=b+m;
        else
            returns(j,i)=2000;
        end
    end

    tmp=wf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    a=std(tmp(50:450,:));
    for j=1:200
        [kWF(i,j),m]=max(tmp(2500:3000,j));
        b=find(tmp(2500+m:4500,j)<a(j),1);
        if ~isempty(b)
            returnsWF(j,i)=b+m;
        else
            returnsWF(j,i)=2000;
        end
    end  
      
    subplot(2,4,i)
    plot(returns(:,i),returnsWF(:,i),'o')
    
    axis([0 2000 0 2000])
    xlabel('Transiency Black')
    ylabel('Transiency White')
    title(['ND', int2str(9-i),])
    line([0 2000],[0 2000],'color','k')
end
data=zeros(200,8,2);
data(:,:,1)=returns;
data(:,:,2)=returnsWF;

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_Z7','data');


%Fig.Z8

for i=1:8
    tmp=bf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);    
    tmp1=wf{i}(:,onOff<0);
    tmp1=tmp1-repmat(mean(tmp1(50:450,:)),4500,1);
    for j=1:200
        cors(i,j)=corr(tmp(501:2500,j),tmp1(2501:4500,j));
        corsWF(i,j)=corr(tmp(2501:4500,j),tmp1(501:2500,j));
    end
end
figure
bar([sum(cors'<-0.2)/2; sum(corsWF'<-0.2)/2]')
xlabel('ND')
ylabel('% of all OFF cells')
title('% of cell with correlation of OFF / ON responses <-0.2')
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
axis([0 9 0 14])
legend('OFF','ON')


%Fig.Z8_alt
for i=1:8
    tmp=bf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);    
    tmp1=wf{i}(:,onOff<0);
    tmp1=tmp1-repmat(mean(tmp1(50:450,:)),4500,1);
    for j=1:200
        cors(i,j)=corr(tmp(501:2500,j),tmp1(2501:4500,j));
        corsON(i,j)=corr(tmp(2501:4500,j),tmp1(501:2500,j));
    end
end
figure
data=zeros(10,4,8);
for i=1:8
    subplot(2,4,i)
    [k,kind]=hist([cors(i,:),-1,1],10);
    k([1,end])=k([1,end])-1;
    k(k<0)=0;
    plot(kind,k/size(cors,2)*100,'-*')
    data(:,1,i)=kind;
    data(:,2,i)=k/size(cors,2)*100;
    
    hold on
    [k,kind]=hist([corsON(i,:),-1,1],10);
    k([1,end])=k([1,end])-1;
    k(k<0)=0;
    plot(kind,k/size(corsON,2)*100,'-*r')
    data(:,3,i)=kind;
    data(:,4,i)=k/size(corsON,2)*100;
    
    title(['ND', int2str(9-i)])
    xlabel('correlation coefficient')
    ylabel('%')
    if i==1
        legend('OFF response','ON response','location','best')
    end
    axis([-1 1 0 80])
end
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_Z8_alt','data');

% calculate white ON response
clear
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all')
nds='87654321';
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/OFFcellsClust_full.txt')
rebSt=0;onSt=0;
on=[0 0 0];reb=[0 0 0];
rebStIs=0;onStIs=0;
rebStWhite=0;onStWhite=0;
onWhite=[0 0 0];rebWhite=[0 0 0];
rebStIsWhite=0;onStIsWhite=0;
clear codes
feat=zeros(12,520);
for i=1:200
    b=OFFcellsClust_full(i,1);
    a=int2str(OFFcellsClust_full(i,2));
    while length(a)<6
        a=['0' a];        
    end
    if length(unique(a(4:6)))==1
        rebSt=rebSt+1;
        if unique(a(4:6))=='1'
            rebStIs=rebStIs+1;
        end
    end

    
    if length(unique(a(1:3)))==1
        onSt=onSt+1;
        if unique(a(1:3))=='1'
            onStIs=onStIs+1;
        end
    end
    
    if a(1)=='1'
        on(1)=on(1)+1;
        feat(1,b)=1;
    end
    if a(2)=='1'
        on(2)=on(2)+1;
        feat(2,b)=1;
    end
    if a(3)=='1'
        on(3)=on(3)+1;
        feat(3,b)=1;
    end
    if a(4)=='1'
        reb(1)=reb(1)+1;
        feat(4,b)=1;
    end
    if a(5)=='1'
        reb(2)=reb(2)+1;
        feat(5,b)=1;
    end
    if a(6)=='1'
        reb(3)=reb(3)+1;
        feat(6,b)=1;
    end    
    aWhite=a;
    
    a=int2str(OFFcellsClust_full(i,3));
    while length(a)<6
        a=['0' a];        
    end
    if length(unique(a(4:6)))==1
        rebStWhite=rebStWhite+1;
        if unique(a(4:6))=='1'
            rebStIsWhite=rebStIsWhite+1;
        end
    end

    
    if length(unique(a(1:3)))==1
        onStWhite=onStWhite+1;
        if unique(a(1:3))=='1'
            onStIsWhite=onStIsWhite+1;
        end
    end
    
    if a(1)=='1'
        onWhite(1)=onWhite(1)+1;
        feat(7,b)=1;
    end
    if a(2)=='1'
        onWhite(2)=onWhite(2)+1;
        feat(8,b)=1;
    end
    if a(3)=='1'
        onWhite(3)=onWhite(3)+1;
        feat(9,b)=1;
    end
    if a(4)=='1'
        rebWhite(1)=rebWhite(1)+1;
        feat(10,b)=1;
    end
    if a(5)=='1'
        rebWhite(2)=rebWhite(2)+1;
        feat(11,b)=1;
    end
    if a(6)=='1'
        rebWhite(3)=rebWhite(3)+1;
        feat(12,b)=1;
    end
    if strcmp(aWhite,a) %if symmetrical, then 1
        ifSym(i)=1;
    else
        ifSym(i)=0;
    end
    codes(i,1:6)=a;
    codes(i,7:12)=aWhite;
end

clear a b c d
for i=0:2
    a(i+1)=sum(~feat(1+i,onOff<0)&~feat(4+i,onOff<0)&~feat(7+i,onOff<0)&~feat(10+i,onOff<0))/2;
    tmp=sum(~feat(1+i,onOff<0)&~feat(7+i,onOff<0))
    b(i+1)=sum(~feat(1+i,onOff<0)&~feat(4+i,onOff<0)&~feat(7+i,onOff<0)&~feat(10+i,onOff<0))/tmp*100;
    c(i+1)=sum(~feat(1+i,onOff<0)&~feat(4+i,onOff<0)&~feat(7+i,onOff<0))/tmp*100;
    d(i+1)=sum(~feat(1+i,onOff<0)&~feat(7+i,onOff<0)&~feat(10+i,onOff<0))/tmp*100;
end

% Fig Z4_white
figure
bar(onWhite/2)
hold on
bar(6:8,rebWhite/2,'facecolor',[0.5 0.1 0])
set(gca,'xtick',[1:3 6:8],'xticklabel',{'ND6','ND5','ND4','ND6','ND5','ND4'})
text(1,-5,'Early ON response')
text(6,-5,'Delayed ON response')
ylabel('% of OFF cells recorded')
title('Fig.Z4. Fraction of cells with early and delayed ON response: white step')

%Fig Z5_white
figure
a=[rebStIsWhite rebStWhite-rebStIsWhite]/2
b=[onStIsWhite onStWhite-onStIsWhite]/2
bar([b; a],'BarLayout','stacked')
set(gca,'xtick',1:2,'xticklabel',{'Early ON response','Delayed ON response'})
ylabel('% of OFF cells recorded') 
axis([0 3 0 70])
legend('always present','never present')
title('Fig.Z5. % of OFF cells with stably present or absent ON response: white step')

data=[onWhite/2 rebWhite/2];
data1=[b; a];
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_Z4_and_Z5_white','data','data1');


sum(ifSym)/2
codes(logical(ifSym),:)

%Fig.Z9 - example of cell with asymmetrical ON responses
figure
nds='87654321'
data=zeros(4500,8,2);
for i=1:8
    subplot(2,4,i)
    tmp=wf{i}(:,129);
    hold off
    plot(tmp,'color','k','linewidth',2)
    data(:,i,1)=tmp;
    
    hold on
    tmp=bf{i}(:,129);
    plot(4701:9200,tmp,'color','k','linewidth',2)
    data(:,i,2)=tmp;
    
    axis([1 9200 0 100])
    line([1 4500],[0,0],'color','k')
    line([4700 9200],[0,0],'color','k')
    
    line([1 500],[55,55]+30,'color','k')
    line([500 2500],[60,60]+30,'color','k')
    line([2500 4500],[55,55]+30,'color','k')
    line([500 500],[55,60]+30,'color','k')
    line([2500 2500],[55,60]+30,'color','k')
    
    line([1 500]+4700,[55,55]+30,'color','k')
    line([500 2500]+4700,[50,50]+30,'color','k')
    line([2500 4500]+4700,[55,55]+30,'color','k')
    line([500 500]+4700,[50,55]+30,'color','k')
    line([2500 2500]+4700,[50,55]+30,'color','k')


    title(['ND',nds(i)])
end

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_Z9','data');


%Fig.Z10 - example of cell without On response at any ND
figure
nds='87654321'
data=zeros(4500,8,2);
for i=1:8
    subplot(2,4,i)
    tmp=wf{i}(:,294);
    tmp=tmp-mean(tmp(50:450,:)); 
    hold off
    plot(tmp,'color','k','linewidth',2)
    data(:,i,1)=tmp;
        
    hold on
    tmp=bf{i}(:,294);
    tmp=tmp-mean(tmp(50:450,:)); 
    plot(4701:9200,tmp,'color','k','linewidth',2)
    data(:,i,2)=tmp;
    
    axis([1 9200 0 100])
    line([1 4500],[0,0],'color','k')
    line([4700 9200],[0,0],'color','k')
    
    line([1 500],[55,55]+30,'color','k')
    line([500 2500],[60,60]+30,'color','k')
    line([2500 4500],[55,55]+30,'color','k')
    line([500 500],[55,60]+30,'color','k')
    line([2500 2500],[55,60]+30,'color','k')
    
    line([1 500]+4700,[55,55]+30,'color','k')
    line([500 2500]+4700,[50,50]+30,'color','k')
    line([2500 4500]+4700,[55,55]+30,'color','k')
    line([500 500]+4700,[50,55]+30,'color','k')
    line([2500 2500]+4700,[50,55]+30,'color','k')

    title(['ND',nds(i)])
end
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_Z10','data');



%Fig.Z15 - example of cell with consistent response change for black and white step
figure
nds='87654321'
data=zeros(4500,8,2);
cnt=37
for i=1:8
    subplot(2,4,i)
    tmp=wf{i}(:,cnt);
    tmp=tmp-mean(tmp(50:450,:));
    hold off
    plot(tmp,'color','k','linewidth',2)
    data(:,i,1)=tmp;
        
    hold on
    tmp=bf{i}(:,cnt);
    tmp=tmp-mean(tmp(50:450,:));
    plot(4701:9200,tmp,'color','k','linewidth',2)
    data(:,i,2)=tmp;
    
    axis([1 9200 0 100])
    line([1 4500],[0,0],'color','k')
    line([4700 9200],[0,0],'color','k')
    
    line([1 500],[55,55]+30,'color','k')
    line([500 2500],[60,60]+30,'color','k')
    line([2500 4500],[55,55]+30,'color','k')
    line([500 500],[55,60]+30,'color','k')
    line([2500 2500],[55,60]+30,'color','k')
    
    line([1 500]+4700,[55,55]+30,'color','k')
    line([500 2500]+4700,[50,50]+30,'color','k')
    line([2500 4500]+4700,[55,55]+30,'color','k')
    line([500 500]+4700,[50,55]+30,'color','k')
    line([2500 2500]+4700,[50,55]+30,'color','k')

    title(['ND',nds(i)])
end
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_Z15','data');


%% Figs.R

% Fig.R
figure
nds='87654321'
data=zeros(4500,2,8);
for i=1:8
    subplot(2,4,i)
    tmp=wf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    hold off
    plot(mean(tmp,2),'color','k','linewidth',2)
    data(:,1,i)=mean(tmp,2);
    
    hold on
    tmp=bf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    plot(4701:9200,mean(tmp,2),'color','k','linewidth',2)
    data(:,2,i)=mean(tmp,2);
    
    axis([1 9200 -20 50])
    line([1 4500],[0,0],'color','k')
    line([4700 9200],[0,0],'color','k')
    
    line([1 500],[40,40],'color','k')
    line([500 2500],[45,45],'color','k')
    line([2500 4500],[40,40],'color','k')
    line([500 500],[40,45],'color','k')
    line([2500 2500],[40,45],'color','k')
    
    line([1 500]+4700,[40,40],'color','k')
    line([500 2500]+4700,[35,35],'color','k')
    line([2500 4500]+4700,[40,40],'color','k')
    line([500 500]+4700,[35,40],'color','k')
    line([2500 2500]+4700,[35,40],'color','k')
    title(['ND',nds(i)])
end
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_R','data');

% R1
figure
data=zeros(4500,3);
for i=3:5
    subplot(1,3,i-2)
    tmp=wf{i}(:,onOff>0);
%     tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    a=std(tmp');
    data(:,i-2)=a;
    plot(a,'k','Linewidth',1)
    title(['ND', nds(i)])
    axis([0 4500 0 40])   
    line([1 500],[32,32],'color','k')
    line([500 2500],[37,37],'color','k')
    line([2500 4500],[32,32],'color','k')
    line([500 500],[32,37],'color','k')
    line([2500 2500],[32,37],'color','k')
    set(gca,'xtick',1000:1000:4000,'xticklabel',{'1','2','3','4'})
    xlabel('time,s')
    ylabel('std')
end
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_R1','data');


% Fig R2

figure
for i=3
    a=wf{i}(:,[190 511 368]);
    plot(a,'linewidth',2)
    line([1 500],[140,140],'color','k')
    line([500 2500],[130,130]+20,'color','k')
    line([2500 4500],[140,140],'color','k')
    line([500 500],[130,140]+10,'color','k')
    line([2500 2500],[140,130]+10,'color','k')
    axis([0 4500 0 160])
    set(gca,'xtick',1000:1000:4000,'xticklabel',{'1','2','3','4'})
    xlabel('time,s')
    ylabel('Firing rate, Hz')
    title('Fig.R2 Example response of 3 ON cells to white step at ND6')
end
data=wf{i}(:,[190 511 368]);
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_R2','data');


% Fig R3
figure
res=cell(8,2);
data=zeros(4500,3,8);
bardata=zeros(3,8);
cellTypes=zeros(3,249);
clear p
for i=1:8
    subplot(2,4,i)
    
    tmp=wf{i}(:,onOff>0);
    varSpont=std(tmp(50:450,:));
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    varSpont(varSpont<1)=1;
    lateSust=mean(tmp(2000:2500,:))>varSpont; % sustained in the last 500ms
    earlyTrans=mean(tmp(800:2500,:))<1&mean(tmp(2000:2500,:))<1; % transient - no response in 500-2000ms
    gaps=mean(tmp(650:1000,:))<mean(tmp(1500:2500,:)); % cells with intermediate inhibition

    trans=earlyTrans;
    inh=gaps&~trans;
    sus=lateSust&~gaps;
    rest=find(~trans&~sus&~inh);
        
    m0=find(trans);
    m=find(inh);
    m1=find(sus);
    
    transMean=mean(tmp(:,m0),2);
    inhMean=mean(tmp(:,m),2);
    susMean=mean(tmp(:,m1),2);
    for k=1:length(rest)
        a=corr(transMean,tmp(:,rest(k)));
        b=corr(inhMean,tmp(:,rest(k)));
        c=corr(susMean,tmp(:,rest(k)));

        [val ind]=max([a b c]);                
        res{i,1}=[res{i,1} val];
        res{i,2}=[res{i,2} ind];
        if val>0.7
            if ind==1
                m0=[m0 rest(k)];
            elseif ind==2
                m=[m rest(k)];
            else
                m1=[m1 rest(k)];
            end
        end            
    end    
    
    transMean=mean(tmp(:,m0),2);
    inhMean=mean(tmp(:,m),2);
    susMean=mean(tmp(:,m1),2);
    data(:,1,i)=transMean;
    data(:,2,i)=inhMean;
    data(:,3,i)=susMean;
    
    if i>2&i<6
        cellTypes(i-2,m0)=1;
        cellTypes(i-2,m)=2;
        cellTypes(i-2,m1)=3;
    end
    
    hold on
    plot(transMean,'color',[0.3 0.8 0.1],'linewidth',2)
    plot(inhMean,'linewidth',2)
    plot(susMean,'r','linewidth',2)
    bar(4650,length(m0)/2.49,180,'facecolor',[0.3 0.8 0.1])
    bar(4850,length(m)/2.49,180,'facecolor','b')
    bar(5050,length(m1)/2.49,180,'facecolor','r')
    
    bardata(1,i)=length(m0)/2.49;
    bardata(2,i)=length(m)/2.49;
    bardata(3,i)=length(m1)/2.49;
    
    
    axis([0 5250 -20 80])
    title(['ND',nds(i)])
    legend(['trans ',int2str(round(length(m0)/2.49)),'%'],['gap ',int2str(round(length(m)/2.49)),'%'],['sust ',int2str(round(length(m1)/2.49)),'%'])
    ylabel('Hz')
    xlabel('time')
    line([1 500],[50,50]-5,'color','k')
    line([500 2500],[55,55]-5,'color','k')
    line([2500 4500],[50,50]-5,'color','k')
    line([500 500],[55,50]-5,'color','k')
    line([2500 2500],[55,50]-5,'color','k')
    (length(m0)+length(m)+length(m1))/2.49
    
    tmp=wf{i}(:,onOff>0);
    tmp=mean(tmp(50:450,:));
    spont(i,1)=mean(tmp(m0));
    spont(i,2)=mean(tmp(m));
    spont(i,3)=mean(tmp(m1));
    spont_std(i,1)=std(tmp(m0))/sqrt(length(m0));
    spont_std(i,2)=std(tmp(m))/sqrt(length(m));
    spont_std(i,3)=std(tmp(m1))/sqrt(length(m1));
    
    p(i,1)=ranksum(tmp(m0),tmp(m));
    p(i,2)=ranksum(tmp(m),tmp(m1));
    p(i,3)=ranksum(tmp(m0),tmp(m1));
end
for i=1:3
    stableCells(i)=sum(sum(cellTypes==i)==3)/2.49;
end
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_R3_R4_R5_R6','data','bardata','spont','spont_std','stableCells','p');

%Fig R4
figure
for i=1:8
    subplot(2,4,i)
    errorbar(100:100:300,spont(i,:),spont_std(i,:),'.k')
    hold on
    bar(100,spont(i,1),90,'facecolor',[0.3 0.8 0.1])
    bar(200,spont(i,2),90,'facecolor','b')
    bar(300,spont(i,3),90,'facecolor','r')
    set(gca,'xtick',100:100:300,'xticklabel',{'trans','gap','sust'})
    ylabel('Mean spontaneous firing rate')
    title(['ND',int2str(9-i)])
    axis([0 400 0 30])
end




%Fig R5
figure
col='bgr'
clear p
subtrspont=0;
for j=1:3
    for i=3:5
        subplot(1,3,j)
        tmp=wf{i}(:,onOff>0);
        varSpont=std(tmp(50:450,:));
        tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
        varSpont(varSpont<1)=1;
        lateSust=mean(tmp(2000:2500,:))>varSpont; % sustained in the last 500ms
        earlyTrans=mean(tmp(800:2500,:))<1&mean(tmp(2000:2500,:))<1; % transient - no response in 500-2000ms
        gaps=mean(tmp(650:1000,:))<mean(tmp(1500:2500,:)); % cells with intermediate inhibition
        
        trans=earlyTrans;
        inh=gaps&~trans;
        sus=lateSust&~gaps;
        rest=find(~trans&~sus&~inh);
        
        m0=find(trans);
        m=find(inh);
        m1=find(sus);
        
        transMean=mean(tmp(:,m0),2);
        inhMean=mean(tmp(:,m),2);
        susMean=mean(tmp(:,m1),2);
        for k=1:length(rest)
            a=corr(transMean,tmp(:,rest(k)));
            b=corr(inhMean,tmp(:,rest(k)));
            c=corr(susMean,tmp(:,rest(k)));
            
            [val ind]=max([a b c]);
            res{i,1}=[res{i,1} val];
            res{i,2}=[res{i,2} ind];
            if val>0.7
                if ind==1
                    m0=[m0 rest(k)];
                elseif ind==2
                    m=[m rest(k)];
                else
                    m1=[m1 rest(k)];
                end
            end
        end
        
       
        
        if ~subtrspont
            tmp=wf{i}(:,onOff>0);
        end
        switch j
            case 1
                tit='transient cells';
                toPlot=mean(tmp(:,m0),2);   
                p(i-2)=length(m0)/2.49;
            case 2
                tit='gap cells';
                toPlot=mean(tmp(:,m),2);   
                p(i-2)=length(m)/2.49;
            case 3
                tit='sustained cells';
                toPlot=mean(tmp(:,m1),2);   
                p(i-2)=length(m1)/2.49;
        end
        
        plot(toPlot,col(i-2),'Linewidth',2);
        hold on
        axis([0 5000 -20 60])
        title(tit)
        
    end
    legend('ND6','ND5','ND4')
    for i=1:3
        bar(4500+i*130,p(i),120,'facecolor',col(i))
    end
end



% Fig.R6
cellTypes=zeros(3,249);
for i=3:5
    
    tmp=wf{i}(:,onOff>0);
    varSpont=std(tmp(50:450,:));
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    varSpont(varSpont<1)=1;
    lateSust=mean(tmp(2000:2500,:))>varSpont; % sustained in the last 500ms
    earlyTrans=mean(tmp(800:2500,:))<1&mean(tmp(2000:2500,:))<1; % transient - no response in 500-2000ms
    gaps=mean(tmp(650:1000,:))<mean(tmp(1500:2500,:)); % cells with intermediate inhibition

    trans=earlyTrans;
    inh=gaps&~trans;
    sus=lateSust&~gaps;
    rest=find(~trans&~sus&~inh);
        
    m0=find(trans);
    m=find(inh);
    m1=find(sus);
    
    transMean=mean(tmp(:,m0),2);
    inhMean=mean(tmp(:,m),2);
    susMean=mean(tmp(:,m1),2);
    for k=1:length(rest)
        a=corr(transMean,tmp(:,rest(k)));
        b=corr(inhMean,tmp(:,rest(k)));
        c=corr(susMean,tmp(:,rest(k)));

        [val ind]=max([a b c]);                
        res{i,1}=[res{i,1} val];
        res{i,2}=[res{i,2} ind];
        if val>0.7
            if ind==1
                m0=[m0 rest(k)];
            elseif ind==2
                m=[m rest(k)];
            else
                m1=[m1 rest(k)];
            end
        end            
    end    
    cellTypes(i-2,m0)=1;
    cellTypes(i-2,m)=2;
    cellTypes(i-2,m1)=3;
end
for i=1:3
    stableCells(i)=sum(sum(cellTypes==i)==3)/2.49;
end
figure
bar(stableCells)
set(gca,'xticklabel',{'trans','gap','sust'})
ylabel('% of all ON cells')
title('ON cells not changing response type')



% Fig.R7
redCT=cellTypes;
i=1;c=1;types=cell(26,1);
while i<size(redCT,2)
    while sum(redCT(:,i)==0)>1
        i=i+1;
    end
    a=redCT-repmat(redCT(:,i),1,size(redCT,2));
    b=sum(a==0)==3;
    types{c}=[find(b) i];
    redCT(:,[find(b) i])=0;
    c=c+1;
    while sum(redCT(:,i))==0&i<size(redCT,2)
        i=i+1;
    end
end
redCT=cellTypes;
clear k
c=1;
for i=1:26
    if length(types{i})>10
        k(1,c)=length(types{i});
        k(2:4,c)=redCT(:,types{i}(1))';
        k(5,c)=i;
        c=c+1;
    end
end

for i=1:7
    tmp=wf{3}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    subplot(3,7,i)
    plot(mean(tmp(:,types{k(5,i)}),2))
end
for i=1:7
    tmp=wf{4}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    subplot(3,7,i+7)
    plot(mean(tmp(:,types{k(5,i)}),2))
end
for i=1:7
    tmp=wf{5}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    subplot(3,7,i+14)
    plot(mean(tmp(:,types{k(5,i)}),2))
end

% amount of cells keeping elevated late response at all 3 light levels
m=0;
for i=1:26
    if sum(cellTypes(:,types{i}(1))==3)+sum(cellTypes(:,types{i}(1))==2)==3
        m=m+length(types{i});
    end
end
m=m/2.49;% in %




% fig R8 (like Z6_alt)
clear p h n sust_fav
data=cell(8,1);
figure
for i=1:8
    
    tmp=wf{i}(:,onOff>0);
    m=sum(tmp(500:1000,:))./sum(tmp(500:2500,:));    
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    [a lat]=max(tmp(500:1000,:));
    
    tmp=bf{i}(:,onOff>0);
    mOn=sum(tmp(2500:3000,:))./sum(tmp(2500:4500,:));
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    [aOn latOn]=max(tmp(2500:3000,:));
      
    %amplitude
    subplot(2,4,i)
    t=aOn./a;
    k=latOn-lat;
    b=isinf(t);
    plot(k(~b),t(~b),'o')
    data{i}=[k(~b);t(~b)];
    
    axis([-50 50 -0 2])
    xlabel('latency delta,ms. Black-White')
    ylabel('amplitude ratio. Black/White')
    title(['ND', int2str(9-i),])
    line([0,0],[0,2],'color','k')
    line([-50 50],[1,1],'color','k')

    m=m(~isnan(m));
    mOn=mOn(~isnan(mOn));
    p(i)=ranksum(a,aOn);
    h(i)=ranksum(lat,latOn);
    n(i)=ranksum(m,mOn);
 
end

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_R8','data')

% fig R9 (like fig Z7)
figure
for i=1:8
    
    tmp=wf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    a=std(tmp(50:450,:));
    for j=1:249
        [k(i,j),m]=max(tmp(500:1000,j));
        b=find(tmp(500+m:2500,j)<a(j),1);
        if ~isempty(b)
            returns(j,i)=b+m;
        else
            returns(j,i)=2000;
        end
    end

    tmp=bf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    a=std(tmp(50:450,:));
    for j=1:249
        [kWF(i,j),m]=max(tmp(2500:3000,j));
        b=find(tmp(2500+m:4500,j)<a(j),1);
        if ~isempty(b)
            returnsBF(j,i)=b+m;
        else
            returnsBF(j,i)=2000;
        end
    end  
      
    subplot(2,4,i)
    plot(returns(:,i),returnsBF(:,i),'o')

    axis([0 2000 0 2000])
    xlabel('Transiency White')
    ylabel('Transiency Black')
    title(['ND', int2str(9-i),])
    line([0 2000],[0 2000],'color','k')
end
data=zeros(249,8,2);
data(:,:,1)=returns;
data(:,:,2)=returnsBF;

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_R9','data');


% R10 like R9 for OFF response
figure
clear returns returnsBF
for i=1:8
    
    tmp=wf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    a=std(tmp(50:450,:));
    for j=1:249
        [k(i,j),m]=min(tmp(2500:3000,j));
        b=find(tmp(2500+m:4500,j)>-a(j),1);
        if ~isempty(b)
            returns(j,i)=b+m;
        else
            returns(j,i)=2000;
        end
    end

    tmp=bf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    a=std(tmp(50:450,:));
    for j=1:249
        [kBF(i,j),m]=min(tmp(500:1000,j));
        b=find(tmp(500+m:2500,j)>-a(j),1);
        if ~isempty(b)
            returnsBF(j,i)=b+m;
        else
            returnsBF(j,i)=2000;
        end
    end  
      
    subplot(2,4,i)
    plot(returns(:,i),returnsBF(:,i),'o')

    axis([0 2000 0 2000])
    xlabel('Recovery White')
    ylabel('Recovery Black')
    title(['ND', int2str(9-i),])
    line([0 2000],[0 2000],'color','k')
end
data=zeros(249,8,2);
data(:,:,1)=returns;
data(:,:,2)=returnsBF;

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_R10','data');



% Fig N7

figure
res=zeros(2,8,2);

bord=1500;
for i=1:8
    subplot(2,4,i)
    
    tmp=wf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    a=sum(tmp(500:bord,:))/(bord-500);
%     a=max(tmp(500:bord,:));
    
    tmp=bf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    b=abs(sum(tmp(500:bord,:))/(bord-500)); 
%     b=abs(min(tmp(500:bord,:))); 
    
    t=frMean_HC(onOff>0,i)-frMean_spont(onOff>0,i); 
    
    hold on
    m=a(t>0)-b(t>0);
    res(1,i,1)=nanmean(m);
    res(1,i,2)=nanstd(m)/sqrt(sum(~isnan(m)));
    
    m=a(t<0)-b(t<0);
    res(2,i,1)=nanmean(m);
    res(2,i,2)=nanstd(m)/sqrt(sum(~isnan(m)));
    
    plot(mean(a(t>0)),mean(b(t>0)),'o')
    plot(mean(a(t<0)),mean(b(t<0)),'or')
    line([0 50],[0 50],'color','k')
    xlabel('exc')
    ylabel('inh')

end
figure
for i=1:8
    subplot(2,4,i)
    errorbar(res(:,i,1),res(:,i,2),'.k')
    hold on
    bar(res(:,i,1))
end
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_N7','res');


% Fig.R15
res=zeros(2000,8,2);
for i=1:8
    subplot(2,4,i)
    t=frMean_HC(onOff>0,i)-frMean_spont(onOff>0,i);
    
    tmp=wf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    tmp1=bf{i}(:,onOff>0);
    tmp1=tmp1-repmat(mean(tmp1(50:450,:)),4500,1);
    
    b=abs(sum(tmp(500:bord,:))/(bord-500));
    for bord=501:2500
        a=sum(tmp(500:bord,:))/(bord-500);
        b=abs(sum(tmp1(500:bord,:))/(bord-500));
        m=a(t>0)-b(t>0);
        res(bord-500,i,1)=nanmean(m);
        m=a(t<0)-b(t<0);
        res(bord-500,i,2)=nanmean(m);
    end
    
end
% save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_R15','res')
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_R15','res')
figure
for i=1:8
    subplot(2,4,i)
    plot(res(:,i,1))
    hold on
    plot(res(:,i,2),'r')
    line([0,2000],[0 0],'color','k')
end



% Fig R11
figure
res=cell(8,2);
cellTypes=zeros(3,249);
data=zeros(4500,3,8);
bardata=zeros(3,8);
for i=1:8
    subplot(2,4,i)
    
    tmp=bf{i}(:,onOff>0);
    varSpont=std(tmp(50:450,:));
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    varSpont(varSpont<1)=1;
    lateSust=mean(tmp(4000:4500,:))>varSpont; % sustained in the last 500ms
    earlyTrans=mean(tmp(2800:4500,:))<1&mean(tmp(2000:2500,:))<1; % transient - no response in 500-2000ms
    gaps=mean(tmp(2650:3000,:))<mean(tmp(3500:4500,:)); % cells with intermediate inhibition

    trans=earlyTrans;
    inh=gaps&~trans;
    sus=lateSust&~gaps;
    rest=find(~trans&~sus&~inh);
        
    m0=find(trans);
    m=find(inh);
    m1=find(sus);
    
    transMean=mean(tmp(:,m0),2);
    inhMean=mean(tmp(:,m),2);
    susMean=mean(tmp(:,m1),2);
    for k=1:length(rest)
        a=corr(transMean,tmp(:,rest(k)));
        b=corr(inhMean,tmp(:,rest(k)));
        c=corr(susMean,tmp(:,rest(k)));

        [val ind]=max([a b c]);                
        res{i,1}=[res{i,1} val];
        res{i,2}=[res{i,2} ind];
        if val>0.7
            if ind==1
                m0=[m0 rest(k)];
            elseif ind==2
                m=[m rest(k)];
            else
                m1=[m1 rest(k)];
            end
        end            
    end    
    
    transMean=mean(tmp(:,m0),2);
    inhMean=mean(tmp(:,m),2);
    susMean=mean(tmp(:,m1),2);
    
    data(:,1,i)=transMean;
    data(:,2,i)=inhMean;
    data(:,3,i)=susMean;
    
    if i>2&i<6
        cellTypes(i-2,m0)=1;
        cellTypes(i-2,m)=2;
        cellTypes(i-2,m1)=3;
    end
    
    
    hold on
    plot(transMean,'color',[0.3 0.8 0.1],'linewidth',2)
    plot(inhMean,'linewidth',2)
    plot(susMean,'r','linewidth',2)
    bar(4650,length(m0)/2.49,180,'facecolor',[0.3 0.8 0.1])
    bar(4850,length(m)/2.49,180,'facecolor','b')
    bar(5050,length(m1)/2.49,180,'facecolor','r')
    
    bardata(1,i)=length(m0)/2.49;
    bardata(2,i)=length(m)/2.49;
    bardata(3,i)=length(m1)/2.49;
    
    axis([0 5250 -20 80])
    title(['ND',nds(i)])
    legend(['trans ',int2str(round(length(m0)/2.49)),'%'],['gap ',int2str(round(length(m)/2.49)),'%'],['sust ',int2str(round(length(m1)/2.49)),'%'])
    ylabel('Hz')
    xlabel('time')
    line([1 500],[50,50]-5,'color','k')
    line([500 2500],[55,55]-5,'color','k')
    line([2500 4500],[50,50]-5,'color','k')
    line([500 500],[55,50]-5,'color','k')
    line([2500 2500],[55,50]-5,'color','k')
    (length(m0)+length(m)+length(m1))/2.49
    
    tmp=bf{i}(:,onOff>0);
    tmp=mean(tmp(50:450,:));
    spont(i,1)=mean(tmp(m0));
    spont(i,2)=mean(tmp(m));
    spont(i,3)=mean(tmp(m1));
    spont_std(i,1)=std(tmp(m0))/sqrt(length(m0));
    spont_std(i,2)=std(tmp(m))/sqrt(length(m));
    spont_std(i,3)=std(tmp(m1))/sqrt(length(m1));
    p(i,1)=ranksum(tmp(m0),tmp(m));
    p(i,2)=ranksum(tmp(m),tmp(m1));
    p(i,3)=ranksum(tmp(m0),tmp(m1));
end
for i=1:3
    stableCells(i)=sum(sum(cellTypes==i)==3)/2.49;
end
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_R11','data','bardata','spont','spont_std','stableCells');




% Fig.R12 (like Fig.Z8_alt)
for i=1:8
    tmp=bf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);    
    tmp1=wf{i}(:,onOff>0);
    tmp1=tmp1-repmat(mean(tmp1(50:450,:)),4500,1);
    for j=1:200
        cors(i,j)=corr(tmp(501:2500,j),tmp1(2501:4500,j)); % OFF response
        corsON(i,j)=corr(tmp(2501:4500,j),tmp1(501:2500,j));
    end
end

figure
data=zeros(10,4,8);
for i=1:8
    subplot(2,4,i)
    [k,kind]=hist([cors(i,:),-1,1],10);
    k([1,end])=k([1,end])-1;
    k(k<0)=0;
    plot(kind,k/size(cors,2)*100,'-*')
    data(:,1,i)=kind;
    data(:,2,i)=k/size(cors,2)*100;
    
    hold on
    [k,kind]=hist([corsON(i,:),-1,1],10);
    k([1,end])=k([1,end])-1;
    k(k<0)=0;
    plot(kind,k/size(corsON,2)*100,'-*r')
    data(:,3,i)=kind;
    data(:,4,i)=k/size(corsON,2)*100;
    
    title(['ND', int2str(9-i)])
    xlabel('correlation coefficient')
    ylabel('%')
    if i==1
        legend('OFF response','ON response','location','best')
    end
    axis([-1 1 0 80])
end
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_R12','data');






% Fig.. like R3 for OFF cells
figure
res=cell(8,2);
cellTypes=zeros(3,200);
for i=1:8
    subplot(2,4,i)
    
    tmp=bf{i}(:,onOff<0);
    varSpont=std(tmp(50:450,:));
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    varSpont(varSpont<1)=1;
    lateSust=mean(tmp(2000:2500,:))>varSpont; % sustained in the last 500ms
    earlyTrans=mean(tmp(800:2500,:))<1&mean(tmp(2000:2500,:))<1; % transient - no response in 500-2000ms
    gaps=mean(tmp(650:1000,:))<mean(tmp(1500:2500,:)); % cells with intermediate inhibition

    trans=earlyTrans;
    inh=gaps&~trans;
    sus=lateSust&~gaps;
    rest=find(~trans&~sus&~inh);
        
    m0=find(trans);
    m=find(inh);
    m1=find(sus);
    
    transMean=mean(tmp(:,m0),2);
    inhMean=mean(tmp(:,m),2);
    susMean=mean(tmp(:,m1),2);
    for k=1:length(rest)
        a=corr(transMean,tmp(:,rest(k)));
        b=corr(inhMean,tmp(:,rest(k)));
        c=corr(susMean,tmp(:,rest(k)));

        [val ind]=max([a b c]);                
        res{i,1}=[res{i,1} val];
        res{i,2}=[res{i,2} ind];
        if val>0.7
            if ind==1
                m0=[m0 rest(k)];
            elseif ind==2
                m=[m rest(k)];
            else
                m1=[m1 rest(k)];
            end
        end            
    end    
    
    transMean=mean(tmp(:,m0),2);
    inhMean=mean(tmp(:,m),2);
    susMean=mean(tmp(:,m1),2);
    if i>2&i<6
        cellTypes(i-2,m0)=1;
        cellTypes(i-2,m)=2;
        cellTypes(i-2,m1)=3;
    end
    
    hold on
    plot(transMean,'color',[0.3 0.8 0.1],'linewidth',2)
    plot(inhMean,'linewidth',2)
    plot(susMean,'r','linewidth',2)
    bar(4650,length(m0)/2.49,180,'facecolor',[0.3 0.8 0.1])
    bar(4850,length(m)/2.49,180,'facecolor','b')
    bar(5050,length(m1)/2.49,180,'facecolor','r')
    
    axis([0 5250 -20 80])
    title(['ND',nds(i)])
    legend(['trans ',int2str(round(length(m0)/2)),'%'],['gap ',int2str(round(length(m)/2)),'%'],['sust ',int2str(round(length(m1)/2)),'%'])
    ylabel('Hz')
    xlabel('time')
    line([1 500],[50,50]-5,'color','k')
    line([500 2500],[55,55]-5,'color','k')
    line([2500 4500],[50,50]-5,'color','k')
    line([500 500],[55,50]-5,'color','k')
    line([2500 2500],[55,50]-5,'color','k')
    (length(m0)+length(m)+length(m1))/2
    
    tmp=bf{i}(:,onOff<0);
    tmp=mean(tmp(50:450,:));
    spont(i,1)=mean(tmp(m0));
    spont(i,2)=mean(tmp(m));
    spont(i,3)=mean(tmp(m1));
    spont_std(i,1)=std(tmp(m0))/sqrt(length(m0));
    spont_std(i,2)=std(tmp(m))/sqrt(length(m));
    spont_std(i,3)=std(tmp(m1))/sqrt(length(m1));
    p(i,1)=ranksum(tmp(m0),tmp(m));
    p(i,2)=ranksum(tmp(m),tmp(m1));
    p(i,3)=ranksum(tmp(m0),tmp(m1));
end

for i=1:3
    stableCells(i)=sum(sum(cellTypes==i)==3)/2;
end

%% Fig P

figure
nds='87654321'
data=zeros(4500,4,8);
for i=1:8
    subplot(2,4,i)
    tmp=wf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    plot(mean(tmp,2),'color','r','linewidth',2)
    hold on
    data(:,1,i)=mean(tmp,2);
    
    tmp=wf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    plot(mean(tmp,2),'color','b','linewidth',2)
    data(:,2,i)=mean(tmp,2);
    
    if i==1
        legend('ON','OFF')
    end
    tmp=bf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    plot(4701:9200,mean(tmp,2),'color','r','linewidth',2)
    data(:,3,i)=mean(tmp,2);
    
    tmp=bf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    plot(4701:9200,mean(tmp,2),'color','b','linewidth',2)
    data(:,4,i)=mean(tmp,2);

    axis([1 9200 -20 65])
    line([1 4500],[0,0],'color','k')
    line([4700 9200],[0,0],'color','k')
    
    line([1 500],[55,55],'color','k')
    line([500 2500],[60,60],'color','k')
    line([2500 4500],[55,55],'color','k')
    line([500 500],[55,60],'color','k')
    line([2500 2500],[55,60],'color','k')
    
    line([1 500]+4700,[55,55],'color','k')
    line([500 2500]+4700,[50,50],'color','k')
    line([2500 4500]+4700,[55,55],'color','k')
    line([500 500]+4700,[50,55],'color','k')
    line([2500 2500]+4700,[50,55],'color','k')
    title(['ND',nds(i)])    
end
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_P','data');


%% Fig.N3
cd('/mnt/muench_data/user/alexandra/scripts')
clear
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all')
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summaryzc_late')
boundaries=[155 130 115 85 75 75 85 85];
figure
tr=1.5
for i=3:5
    pol=onOff>0;
    a=zc_HC(i,:)-zc_LC(i,:);
    cond=zc_LC(i,:)>0&zc_HC(i,:)>0&bothContr&zc_LC(i,:)~=boundaries(i)&zc_HC(i,:)~=boundaries(i)&a<-tr;
    a=reshape(HC(:,i,pol&cond),500,sum(pol&cond));
    a=nanmean(a./repmat(sum(abs(a)),500,1),2);
    a=a./sum(abs(a));
    hold on
    plot(a/max(a),'linewidth',2)

    a=reshape(LC(:,i,pol&cond),500,sum(pol&cond));
    a=nanmean(a./repmat(sum(abs(a)),500,1),2);
    a=a./sum(abs(a));
    plot(a/max(a),'r','linewidth',2)
end


% figure
for i=3:5
    pol=onOff<0;
    a=zc_HC(i,:)-zc_LC(i,:);
    cond=zc_LC(i,:)>0&zc_HC(i,:)>0&bothContr&zc_LC(i,:)~=boundaries(i)&zc_HC(i,:)~=boundaries(i)&a<-tr;
    a=reshape(HC(:,i,pol&cond),500,sum(pol&cond));
    a=nanmean(a./repmat(sum(abs(a)),500,1),2);
    a=a./sum(abs(a));
    hold on
    plot(a/min(a),'c','linewidth',2)
    
    a=reshape(LC(:,i,pol&cond),500,sum(pol&cond));
    a=nanmean(a./repmat(sum(abs(a)),500,1),2);
    a=a./sum(abs(a));
    plot(a/min(a),'m','linewidth',2)
end



figure
tr=1.5
for i=3:5
    pol=onOff>0;
    a=zc_HC(i,:)-zc_LC(i,:);
    cond=zc_LC(i,:)>0&zc_HC(i,:)>0&bothContr&zc_LC(i,:)~=boundaries(i)&zc_HC(i,:)~=boundaries(i)&a<-tr;
    a=reshape(HC(:,i,pol&cond),500,sum(pol&cond));
    a=nanmean(a./repmat(sum(abs(a)),500,1),2);
    a=a./sum(abs(a));
    hold on
    plot(a/max(a),'linewidth',2)

    a=reshape(LC(:,i,pol&cond),500,sum(pol&cond));
    a=nanmean(a./repmat(sum(abs(a)),500,1),2);
    a=a./sum(abs(a));
    plot(a/max(a),'r','linewidth',2)
end


% figure
for i=3:5
    pol=onOff>0;
    a=zc_HC(i,:)-zc_LC(i,:);
    cond=zc_LC(i,:)>0&zc_HC(i,:)>0&bothContr&zc_LC(i,:)~=boundaries(i)&zc_HC(i,:)~=boundaries(i)&a>tr;
    a=reshape(HC(:,i,pol&cond),500,sum(pol&cond));
    a=nanmean(a./repmat(sum(abs(a)),500,1),2);
    a=a./sum(abs(a));
    hold on
    plot(a/max(a),'c','linewidth',2)

    a=reshape(LC(:,i,pol&cond),500,sum(pol&cond));
    a=nanmean(a./repmat(sum(abs(a)),500,1),2);
    a=a./sum(abs(a));
    plot(a/max(a),'m','linewidth',2)
end








cd('/mnt/muench_data/user/alexandra/scripts')
clear
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all')
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summaryzc_late')
boundaries=[155 130 115 85 75 75 85 85];
ranges=[2 3 3 4 2;3 4 5 5 5];
figure
cc=1;
for polc=1:2
    for j=1:5
        subplot(2,5,cc)
        hold on
        k=ones(1,520);
        for i=[ranges(1,j) ranges(2,j)]
            if polc==1
                pol=onOff<0;
            else
                pol=onOff>0;
            end
            a=zc_HC(i,:)-zc_LC(i,:);
            cond=zc_LC(i,:)>0&zc_HC(i,:)>0&bothContr&zc_LC(i,:)~=boundaries(i)&zc_HC(i,:)~=boundaries(i)&abs(a)<30;
            k(~(pol&cond))=0;
        end
        ll=sum(k);
        k=logical(k);
        a=zc_HC(ranges(1,j),k)-zc_LC(ranges(1,j),k);
        a=[a; zc_HC(ranges(2,j),k)-zc_LC(ranges(2,j),k)];
        b=a(1,:)-a(2,:);
        plot(a(:,b>=0),'color',[1 0.5 0.5]*0.7);
        plot(a(:,b<0),'color',[0.5 0.5 1]*0.7);
        axis([0.8 2.2 -30 20])
        plot(mean(a'),'k','lineWidth',4)
        line([0 3],[0 0],'color','k','lineWidth',2)
        title([int2str(9-ranges(1,j)),'  ',int2str(9-ranges(2,j)),'  ',int2str(polc-2),'   ',int2str(ll), '  ', int2str(sum(b>=0)),'  ', int2str(sum(b<0))])
        cc=cc+1;
    end
end



ranges=[2 3 3 4 2;3 4 5 5 5]
figure
cc=1;
for polc=1:2
    for j=1:5
        subplot(2,5,cc)
        hold on
        k=ones(1,520);
        for i=[ranges(1,j) ranges(2,j)]
            if polc==1
                pol=onOff<0;
            else
                pol=onOff>0;
            end
            a=peak_HC(i,:)-peak_LC(i,:);
            cond=peak_LC(i,:)>0&peak_HC(i,:)>0&bothContr&zc_LC(i,:)~=boundaries(i)&zc_HC(i,:)~=boundaries(i)&abs(a)<30;
            k(~(pol&cond))=0;
        end
        ll=sum(k);
        k=logical(k);
        a=peak_HC(ranges(1,j),k)-peak_LC(ranges(1,j),k);
        a=[a; peak_HC(ranges(2,j),k)-peak_LC(ranges(2,j),k)];
        b=a(1,:)-a(2,:);
        plot(a(:,b>=0),'color',[1 0.5 0.5]*0.7);
        plot(a(:,b<0),'color',[0.5 0.5 1]*0.7);
        axis([0.8 2.2 -30 20])

        plot(mean(a'),'k','lineWidth',4)
        line([0 3],[0 0],'color','k','lineWidth',2)
        title([int2str(9-ranges(1,j)),'  ',int2str(9-ranges(2,j)),'  ',int2str(polc-2),'   ',int2str(ll), '  ', int2str(sum(b>=0)),'  ', int2str(sum(b<0))])
        cc=cc+1;
    end
end







ranges=[2 3 3 4 2;3 4 5 5 5]
figure
cc=1;
for polc=1:2
    for j=1:5
        subplot(2,5,cc)
        hold on
        k=ones(1,520);
        for i=[3 4 5]
            if polc==1
                pol=onOff<0;
            else
                pol=onOff>0;
            end
            a=peak_HC(i,:)-peak_LC(i,:);
            cond=peak_LC(i,:)>0&peak_HC(i,:)>0&bothContr&zc_LC(i,:)~=boundaries(i)&zc_HC(i,:)~=boundaries(i)&abs(a)<30;
            k(~(pol&cond))=0;
        end
        ll=sum(k);
        k=logical(k);
        a=peak_HC(3,k)-peak_LC(3,k);
        a=[a; peak_HC(4,k)-peak_LC(4,k)];
        a=[a; peak_HC(5,k)-peak_LC(5,k)];
        figure
        plot(a)

%         plot(a(:,b>=0),'color',[1 0.5 0.5]*0.7);
%         plot(a(:,b<0),'color',[0.5 0.5 1]*0.7);
        axis([0.8 2.2 -30 20])

        plot(mean(a'),'k','lineWidth',4)
        line([0 3],[0 0],'color','k','lineWidth',2)
        title([int2str(9-ranges(1,j)),'  ',int2str(9-ranges(2,j)),'  ',int2str(polc-2),'   ',int2str(ll), '  ', int2str(sum(b>=0)),'  ', int2str(sum(b<0))])
        cc=cc+1;
    end
end






figure
cc=1;
for polc=1:2
    for j=1:5
        subplot(2,5,cc)
        hold on
        k=ones(1,520);
        for i=[ranges(1,j) ranges(2,j)]
            if polc==1
                pol=onOff<0;
            else
                pol=onOff>0;
            end
            a=zc_HC(i,:)-zc_LC(i,:);
            cond=zc_LC(i,:)>0&zc_HC(i,:)>0&bothContr&zc_LC(i,:)~=boundaries(i)&zc_HC(i,:)~=boundaries(i)&abs(a)<30;
            k(~(pol&cond))=0;
        end
        ll=sum(k);
        k=logical(k);
        a=zc_HC(ranges(1,j),k)-zc_LC(ranges(1,j),k);
        a=[a; zc_HC(ranges(2,j),k)-zc_LC(ranges(2,j),k)];
        hist([a(1,:)-a(2,:) -50 50],50)
        axis([-40 40 0 15])
%         plot(mean(a'),'k','lineWidth',4)
        line([0 0],[0 10],'color','r','lineWidth',2)
        title([int2str(9-ranges(1,j)),'  ',int2str(9-ranges(2,j)),'  ',int2str(polc-2),'   ',int2str(ll)])
        cc=cc+1;
    end
end





%% Firing rate. Fig. K

clear
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all')

%Fig.K
figure
for i=1:8
    subplot(2,4,i)
    line([0 250],[0 0],'color','k')
    hold on
    b=frMean_HC(onOff>0,i)-frMean_spont(onOff>0,i);
    [v,m]=sort(frMean_spont(onOff>0,i));
    plot(v,b(m),'or','markersize',6)
    
    b=frMean_HC(onOff<0,i)-frMean_spont(onOff<0,i);
    [v,m]=sort(frMean_spont(onOff<0,i));
    plot(v,b(m),'o','markersize',6)
    axis([0 60 -20 40])
    title(['ND',int2str(9-i)])
    xlabel('spont fr')
    ylabel('HC FR - spont FR')
end

%Fig.K
figure
for i=1:8
    subplot(2,4,i)
    line([0 250],[0 0],'color','k')
    hold on
    b=frMean_HC(onOff>0&bothContr,i)-frMean_LC(onOff>0&bothContr,i);
    [v,m]=sort(frMean_spont(onOff>0&bothContr,i));
    plot(v,b(m),'or','markersize',6)
    
    b=frMean_HC(onOff<0&bothContr,i)-frMean_LC(onOff<0&bothContr,i);
    [v,m]=sort(frMean_spont(onOff<0&bothContr,i));
    plot(v,b(m),'o','markersize',6)
    axis([0 60 -10 30])
    title(['ND',int2str(9-i)])
    xlabel('spont fr')
    ylabel('HC FR - spont FR')
end



% Fig.K1
clear l
for i=1:8
    ons=onOff<0&bothContr>0;
    l(i,1)=sum((frMean_HC(ons,i)-frMean_LC(ons,i))<0)/sum(ons);
    ons=onOff>0&bothContr>0;
    l(i,2)=sum((frMean_HC(ons,i)-frMean_LC(ons,i))<0)/sum(ons);
end

hb=bar(l*100);
legend(hb,'OFF','ON')
axis([0 9 0 60])
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND')
ylabel('%')
title('Fig.K1. Amount of cells decreasing firing rate at high contrast')



%% Fig.O
clear
cd('/mnt/muench_data/user/alexandra/scripts')
date='20130302_1'
codeWord='NatMov'
dataRaster=cell(3,1);
dataRaster_supp=zeros(3,1);
data=zeros(4500,2,3);

mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/'];
hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end

% select heka files
file_list=[];
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,codeWord, 'once'))
        file_list=[file_list i];
    end
end
conv_acc=zeros(3,25000);
units=dir([mainpath, 'units/*.mat']);
cnt=26;
ordCnt=387;
units(cnt).name
load([mainpath,'units/',units(cnt).name]);

t=[];

cnt=0;
for i=37:37+17
    spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms       
    conv_tmp=convolved(spikes,40,25000);
    conv_acc(1,:)=conv_acc(1,:)+conv_tmp(121:end-120)/18;
    t=[t spikes+26000*cnt];    
    cnt=cnt+1;
end
p=subplot('Position',[0.05 0.05 0.60 0.24])

dataRaster{1}=t;
dataRaster_supp(1)=cnt;

rasterplot(t,cnt,26000,p);
title('ND6')
t=[];
cnt=0;
for i=55:55+17
    spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms 
    conv_tmp=convolved(spikes,40,25000);
    conv_acc(2,:)=conv_acc(2,:)+conv_tmp(121:end-120)/18;
    t=[t spikes+26000*cnt];    
    cnt=cnt+1;
end

p=subplot('Position',[0.05 0.38 0.60 0.24])

dataRaster{2}=t;
dataRaster_supp(2)=cnt;

rasterplot(t,cnt,26000,p);
title('ND5')
t=[];
cnt=0;
for i=73:73+17
    spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
    conv_tmp=convolved(spikes,40,25000);
    conv_acc(3,:)=conv_acc(3,:)+conv_tmp(121:end-120)/18;
    t=[t spikes+26000*cnt];    
    cnt=cnt+1;
end

p=subplot('Position',[0.05 0.71 0.60 0.24])

dataRaster{3}=t;
dataRaster_supp(3)=cnt;

rasterplot(t,cnt,26000,p);
title('ND4')

load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all','bf','wf')

p(1)=subplot('Position',[0.70 0.05 0.28 0.24])
a=wf{3}(:,ordCnt);
b=bf{3}(:,ordCnt);
plot(a,'color','k')
hold on
plot(4701:9200,b,'color','k')
title('ND6')

data(:,1,1)=a;
data(:,2,1)=b;

p(2)=subplot('Position',[0.70 0.38 0.28 0.24])
a=wf{4}(:,ordCnt);
b=bf{4}(:,ordCnt);
plot(a,'color','k')
hold on
plot(4701:9200,b,'color','k')
title('ND5')

data(:,1,2)=a;
data(:,2,2)=b;


p(3)=subplot('Position',[0.70 0.71 0.28 0.24])
a=wf{5}(:,ordCnt);
b=bf{5}(:,ordCnt);
plot(a,'color','k')
hold on
plot(4701:9200,b,'color','k')
title('ND4')

data(:,1,3)=a;
data(:,2,3)=b;

for i=1:3
    subplot(p(i))
    set(p(i),'xtick',[500 2500 5300 7300],'xticklabel',{'0.5','2.5','0.5','2.5'})
    xlabel('time,s')
    ylabel('Frequency, Hz')
    axis([0 9200 0 125])
        
    line([1 500],[55,55]+55,'color','k')
    line([500 2500],[60,60]+55,'color','k')
    line([2500 4500],[55,55]+55,'color','k')
    line([500 500],[55,60]+55,'color','k')
    line([2500 2500],[55,60]+55,'color','k')
    
    line([1 500]+4700,[55,55]+55,'color','k')
    line([500 2500]+4700,[50,50]+55,'color','k')
    line([2500 4500]+4700,[55,55]+55,'color','k')
    line([500 500]+4700,[50,55]+55,'color','k')
    line([2500 2500]+4700,[50,55]+55,'color','k')
    
end

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_O','dataRaster','dataRaster_supp','data','conv_acc')

%% Fig.E

load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all')
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_gain','hc_gain','lc_gain')
for i=1:8
    for j=find(bothContr)
        [ind,lags]=xcorr(HC(50:300,i,j),LC(50:300,i,j),'coeff');
        [mc,ic]=max(ind);
        delay(j,i)=lags(ic);
        value(j,i)=mc;
        gainContr(j,i)=hc_gain(j,i)./lc_gain(j,i);
        fr(j,i)=frMean_HC(j,i)-frMean_LC(j,i);
    end
end
figure

subplot(3,3,1)
plot(1,1,'b')
hold on
plot(1,1,'r')
legend('Increase FR','Decrease FR')
plot(1,1,'color',[1 1 1])
plnmr=[4 7 2 5 8 3 6 9];
data=zeros(4500,4,8);
for i=1:8
    subplot(3,3,plnmr(i))
    tmp=bf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    a=fr(onOff>0,i);
    
    tmp1=wf{i}(:,onOff>0);
    tmp1=tmp1-repmat(mean(tmp1(50:450,:)),4500,1);
    plot(mean(tmp1(:,a<0),2),'r');    
    hold on
    plot(mean(tmp1(:,a>0),2),'b');
    data(:,1,i)=mean(tmp1(:,a<0),2);
    data(:,2,i)=mean(tmp1(:,a>0),2);
    
    
    plot(4701:9200,mean(tmp(:,a<0),2),'r');
    hold on
    plot(4701:9200,mean(tmp(:,a>0),2),'b');
    data(:,3,i)=mean(tmp(:,a<0),2);
    data(:,4,i)=mean(tmp(:,a>0),2);
    
    axis([0 9200 -25 65])

    title(['ND',int2str(9-i)])
    xlabel('time')
    ylabel('firing rate, Hz')

    line([1 4500],[0,0],'color','k')
    line([4700 9200],[0,0],'color','k')
    
    line([1 500],[55,55],'color','k')
    line([500 2500],[60,60],'color','k')
    line([2500 4500],[55,55],'color','k')
    line([500 500],[55,60],'color','k')
    line([2500 2500],[55,60],'color','k')
    
    line([1 500]+4700,[55,55],'color','k')
    line([500 2500]+4700,[50,50],'color','k')
    line([2500 4500]+4700,[55,55],'color','k')
    line([500 500]+4700,[50,55],'color','k')
    line([2500 2500]+4700,[50,55],'color','k')
end
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_E','data')

%% Fig. E1

clear
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all')
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_gain','hc_gain','lc_gain')
for i=1:8
    for j=find(bothContr)
        [ind,lags]=xcorr(HC(50:300,i,j),LC(50:300,i,j),'coeff');
        [mc,ic]=max(ind);
        delay(j,i)=lags(ic);
        value(j,i)=mc;
        gainContr(j,i)=hc_gain(j,i)./lc_gain(j,i);
        fr(j,i)=frMean_HC(j,i)-frMean_LC(j,i);
    end
end


figure

subplot(3,3,1)
plot(1,1,'b')
hold on
plot(1,1,'r')
legend('Speed up','Slow down')
plot(1,1,'color',[1 1 1])
plnmr=[4 7 2 5 8 3 6 9];
for i=1:8
    subplot(3,3,plnmr(i))
    tmp=bf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    a=delay(onOff>0,i);
    b=value(onOff>0,i);
    d=a<30&a>-30&b>0.6;
    a=-a(d);
    tmp=tmp(:,d);
    
    tmp1=wf{i}(:,onOff>0);
    tmp1=tmp1-repmat(mean(tmp1(50:450,:)),4500,1);
    tmp1=tmp1(:,d);
    
    plot(mean(tmp1(:,a<0),2),'r');
    hold on
    plot(mean(tmp1(:,a>0),2),'b');
    
    plot(4701:9200,mean(tmp(:,a<0),2),'r');
    hold on
    plot(4701:9200,mean(tmp(:,a>0),2),'b');

    axis([0 9200 -25 65])

    title(['ND',int2str(9-i)])
    xlabel('time')
    ylabel('firing rate, Hz')

    line([1 4500],[0,0],'color','k')
    line([4700 9200],[0,0],'color','k')
    
    line([1 500],[55,55],'color','k')
    line([500 2500],[60,60],'color','k')
    line([2500 4500],[55,55],'color','k')
    line([500 500],[55,60],'color','k')
    line([2500 2500],[55,60],'color','k')
    
    line([1 500]+4700,[55,55],'color','k')
    line([500 2500]+4700,[50,50],'color','k')
    line([2500 4500]+4700,[55,55],'color','k')
    line([500 500]+4700,[50,55],'color','k')
    line([2500 2500]+4700,[50,55],'color','k')
end





%% Fig.E2 execute chirp_acc first

load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all')
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_gain','hc_gain','lc_gain')
for i=1:8
    for j=find(bothContr)
        [ind,lags]=xcorr(HC(50:300,i,j),LC(50:300,i,j),'coeff');
        [mc,ic]=max(ind);
        delay(j,i)=lags(ic);
        value(j,i)=mc;
        gainContr(j,i)=hc_gain(j,i)./lc_gain(j,i);
        fr(j,i)=frMean_HC(j,i)-frMean_LC(j,i);
    end
end

figure
chirp_acc=zeros(23000,2,8);
for i=1:8
    subplot(4,2,i)
    tmp=reshape(mean(chirp(:,(i-1)*4+1:i*4,onOff>0&real_chirp'>0),2),23000,sum(onOff>0&real_chirp'>0));
    tmp=tmp-repmat(mean(tmp(50:2500,:)),23000,1);
   
    a=fr(onOff>0&real_chirp'>0,i);

    plot(mean(tmp(:,a>0),2),'b');
    hold on    
    plot(mean(tmp(:,a<0),2),'r');
    chirp_acc(:,1,i)=mean(tmp(:,a>0),2);
    chirp_acc(:,2,i)=mean(tmp(:,a<0),2);

    line([0 22000],[0 0],'color','k')
    axis([2500 21000 -30 40])

    title(['ND',int2str(9-i)])
    xlabel('time')
    ylabel('firing rate, Hz')

end
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/chirp_for_E2.mat','chirp_acc')

%% Fig.E3
figure
chirp_acc=zeros(23000,8);
for i=1:8
    subplot(4,2,i)
    tmp=reshape(mean(chirp(:,(i-1)*4+1:i*4,onOff<0&real_chirp'>0),2),23000,sum(onOff<0&real_chirp'>0));
    tmp=tmp-repmat(mean(tmp(50:2500,:)),23000,1);

    plot(mean(tmp,2),'b');
    chirp_acc(:,i)=mean(tmp,2);
    line([0 22000],[0 0],'color','k')
    axis([2500 21000 -10 65])

    title(['ND',int2str(9-i)])
    xlabel('time')
    ylabel('firing rate, Hz')

end
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/chirp_for_E3.mat','chirp_acc')


figure
chirp_acc=zeros(23000,2,8);
for i=1:8
    subplot(4,2,i)
    tmp=reshape(mean(chirp(:,(i-1)*4+1:i*4,onOff<0&real_chirp'>0),2),23000,sum(onOff<0&real_chirp'>0));
    tmp=tmp-repmat(mean(tmp(50:2500,:)),23000,1);

    
    chirp_acc(:,1,i)=mean(tmp(:,1:2:end),2);
    chirp_acc(:,2,i)=mean(tmp(:,2:2:end),2);
    plot(chirp_acc(:,1,i),'b');
    hold on
    plot(chirp_acc(:,2,i),'r');
    line([0 22000],[0 0],'color','k')
    axis([2500 21000 -10 65])

    title(['ND',int2str(9-i)])
    xlabel('time')
    ylabel('firing rate, Hz')
    corr(chirp_acc(:,1,i),chirp_acc(:,2,i))

end


%% Fig.E4
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all')
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_gain','hc_gain','lc_gain')
for i=1:8
    for j=find(bothContr)
        [ind,lags]=xcorr(HC(50:300,i,j),LC(50:300,i,j),'coeff');
        [mc,ic]=max(ind);
        delay(j,i)=lags(ic);
        value(j,i)=mc;
        gainContr(j,i)=hc_gain(j,i)./lc_gain(j,i);
        fr(j,i)=frMean_HC(j,i)-frMean_LC(j,i);
    end
end
figure

data2plot=zeros(23000,8,3);
for i=1:8
    subplot(4,2,i)
    tmp=reshape(mean(chirp(:,(i-1)*4+1:i*4,onOff>0&real_chirp'>0),2),23000,sum(onOff>0&real_chirp'>0));
    tmp=tmp-repmat(mean(tmp(50:2500,:)),23000,1);
   
    a=fr(onOff>0&real_chirp'>0,i);
    hold on    
    plot(mean(tmp(:,a<0),2),'r');
    data2plot(:,i,1)=mean(tmp(:,a<0),2);
    data2plot(:,i,2)=mean(tmp(:,a>0),2);
        
        
    tmp=reshape(mean(chirp(:,(i-1)*4+1:i*4,onOff<0&real_chirp'>0),2),23000,sum(onOff<0&real_chirp'>0));
    tmp=tmp-repmat(mean(tmp(50:2500,:)),23000,1);    
    plot(mean(tmp,2),'b');
    data2plot(:,i,3)=mean(tmp,2);
    
    line([0 22000],[0 0],'color','k')
    axis([2500 21000 -30 60])

    title(['ND',int2str(9-i)])
    xlabel('time')
    ylabel('firing rate, Hz')

end





%% Fig.E2 execute chirp_acc first

load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all')
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_gain','hc_gain','lc_gain')
for i=1:8
    for j=find(bothContr)
        [ind,lags]=xcorr(HC(50:300,i,j),LC(50:300,i,j),'coeff');
        [mc,ic]=max(ind);
        delay(j,i)=lags(ic);
        value(j,i)=mc;
        gainContr(j,i)=hc_gain(j,i)./lc_gain(j,i);
        fr(j,i)=frMean_HC(j,i)-frMean_LC(j,i);
    end
end





figure

for i=3:5
    subplot(3,1,i-2)
    
    tmp=bf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    a=std(tmp(50:450,:));
    for j=1:200
        [k(i,j),m]=max(tmp(500:1000,j));
        b=find(tmp(500+m:2500,j)<a(j),1);
        if ~isempty(b)
            returns(j,i)=b+m;
        else
            returns(j,i)=2000;
        end
    end
c1=returns(:,i)'<500;
c2=returns(:,i)'>1950;
%     tmp=reshape(mean(chirp(:,(i-1)*4+1:i*4,c1&onOff<0&real_chirp'>0),2),23000,sum(c1&onOff<0&real_chirp'>0));
%     tmp=tmp-repmat(mean(tmp(50:2500,:)),23000,1);
%     plot(mean(tmp,2),'r');
%     hold on
%     tmp=reshape(mean(chirp(:,(i-1)*4+1:i*4,c2&onOff<0&real_chirp'>0),2),23000,sum(c2&onOff<0&real_chirp'>0));
%     tmp=tmp-repmat(mean(tmp(50:2500,:)),23000,1);
%     plot(mean(tmp,2),'b');
%     
    k=reshape(HC(:,i,onOff<0),500,200);
    k=k./repmat(sum(abs(k)),500,1);
    plot(mean(k(:,c1),2),'r');
    hold on
    plot(mean(k(:,c2),2),'b');
    
%     tmp=bf{i};
%     tmp=tmp-repmat(mean(tmp(50:1950,:)),4500,1);    
% 
%     plot(mean(tmp(:,c1),2),'b');
%     hold on    
%     plot(mean(tmp(:,c2),2),'r');
%     
%     tmp=wf{i};
%     tmp=tmp-repmat(mean(tmp(50:1950,:)),4500,1);
%     
%     plot(4701:9200,mean(tmp(:,feat(i-2,:)==1&onOff<0),2),'b');
%     hold on    
%     plot(4701:9200,mean(tmp(:,feat(i-2,:)==0&onOff<0),2),'r');   
    
%     plot(mean(tmp(:,b<0),2),'g');
% 
%     line([0 22000],[0 0],'color','k')
%     axis([2500 21000 -40 50])

    title(['ND',int2str(9-i)])
    xlabel('time')
    ylabel('firing rate, Hz')

end


figure
for i=1:8
    
    tmp=bf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    a=std(tmp(50:450,:));
    for j=1:200
        [k(i,j),m]=max(tmp(500:1000,j));
        b=find(tmp(500+m:2500,j)<a(j),1);
        if ~isempty(b)
            returns(j,i)=b+m;
        else
            returns(j,i)=2000;
        end
    end

    tmp=wf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    a=std(tmp(50:450,:));
    for j=1:200
        [kWF(i,j),m]=max(tmp(2500:3000,j));
        b=find(tmp(2500+m:4500,j)<a(j),1);
        if ~isempty(b)
            returnsWF(j,i)=b+m;
        else
            returnsWF(j,i)=2000;
        end
    end  
      
    subplot(2,4,i)
    plot(returns(:,i),returnsWF(:,i),'o')

    axis([0 2000 0 2000])
    xlabel('Transiency Black')
    ylabel('Transiency White')
    title(['ND', int2str(9-i),])
    line([0 2000],[0 2000],'color','k')
end



%% Fig X (cnt=176)
cnt=176
data=zeros(4500,2,8);
for i=1:8
    data(:,1,i)=bf{i}(:,cnt);
    data(:,2,i)=wf{i}(:,cnt);
end
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_X','data')

%% Fig Y (cnt=114)
cnt=114
data=zeros(4500,2,8);
for i=1:8
    data(:,1,i)=bf{i}(:,cnt);
    data(:,2,i)=wf{i}(:,cnt);
end
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_Y','data')
 