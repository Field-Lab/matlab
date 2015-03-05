%% 20130220
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130220'
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

FiringRate=zeros(65000,size(SpikeCount,1),size(SpikeCount,2));
for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,65000);
            FiringRate(:,i,cnt)=tmp(121:end-120);
            tmp=FiringRate(:,i,cnt);

            frSTD_HC(i,cnt)=std(tmp(4000:55000));
            frMean_HC(i,cnt)=mean(tmp(4000:55000));

            frSTD_LC(i,cnt)=-1;
            frMean_LC(i,cnt)=-1;
            frspont(i,cnt)=mean(tmp(200:1950));
            frSTDspont(i,cnt)=std(tmp(200:1950));
        end
    end
end

formula=11:9:82;
save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frspont', 'frSTDspont','frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','formula')

%% 20130220_1
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130220_1'
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

FiringRate=zeros(65000,size(SpikeCount,1),size(SpikeCount,2));
for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,65000);
            FiringRate(:,i,cnt)=tmp(121:end-120);
            tmp=FiringRate(:,i,cnt);

            frSTD_HC(i,cnt)=std(tmp(4000:55000));
            frMean_HC(i,cnt)=mean(tmp(4000:55000));

            frSTD_LC(i,cnt)=-1;
            frMean_LC(i,cnt)=-1;
            frspont(i,cnt)=mean(tmp(200:1950));
            frSTDspont(i,cnt)=std(tmp(200:1950));
        end
    end
end

formula=10:9:81;
save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frspont', 'frSTDspont','frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','formula')

%% 20130224
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130224'
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

FiringRate=zeros(65000,size(SpikeCount,1),size(SpikeCount,2));
for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,65000);
            FiringRate(:,i,cnt)=tmp(121:end-120);
            tmp=FiringRate(:,i,cnt);

            frSTD_HC(i,cnt)=std(tmp(4000:55000));
            frMean_HC(i,cnt)=mean(tmp(4000:55000));

            frSTD_LC(i,cnt)=-1;
            frMean_LC(i,cnt)=-1;
            frspont(i,cnt)=mean(tmp(200:1950));
            frSTDspont(i,cnt)=std(tmp(200:1950));
        end
    end
end

formula=10:9:81;
save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frspont', 'frSTDspont','frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','formula')


%% 20130225
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130225'
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

FiringRate=zeros(65000,size(SpikeCount,1),size(SpikeCount,2));
for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,65000);
            FiringRate(:,i,cnt)=tmp(121:end-120);
            tmp=FiringRate(:,i,cnt);

            frSTD_HC(i,cnt)=std(tmp(4000:55000));
            frMean_HC(i,cnt)=mean(tmp(4000:55000));

            frSTD_LC(i,cnt)=-1;
            frMean_LC(i,cnt)=-1;
            frspont(i,cnt)=mean(tmp(200:1950));
            frSTDspont(i,cnt)=std(tmp(200:1950));
        end
    end
end

formula=10:9:81;
save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frspont', 'frSTDspont','frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','formula')

%% 20130226
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130226'
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

FiringRate=zeros(65000,size(SpikeCount,1),size(SpikeCount,2));
for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,65000);
            FiringRate(:,i,cnt)=tmp(121:end-120);
            tmp=FiringRate(:,i,cnt);

            frSTD_HC(i,cnt)=std(tmp(4000:55000));
            frMean_HC(i,cnt)=mean(tmp(4000:55000));

            frSTD_LC(i,cnt)=-1;
            frMean_LC(i,cnt)=-1;
            frspont(i,cnt)=mean(tmp(200:1950));
            frSTDspont(i,cnt)=std(tmp(200:1950));
        end
    end
end

formula=10:9:81;
save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frspont', 'frSTDspont','frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','formula')


%% 20130301
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130301'
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

FiringRate=zeros(65000,size(SpikeCount,1),size(SpikeCount,2));
for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,65000);
            FiringRate(:,i,cnt)=tmp(121:end-120);
            tmp=FiringRate(:,i,cnt);
            for mc=1:6
                take=correctedProtocols(600*(mc-1)+1,1,i)+2000:correctedProtocols(600*mc,1,i)-500;
                tmp1(:,mc)=tmp(take(1:7000));
            end
            tmp=reshape(tmp1(:,1:2:end),21000,1);
            frSTD_HC(i,cnt)=std(tmp);
            frMean_HC(i,cnt)=mean(tmp);
                        
            tmp=reshape(tmp1(:,2:2:end),21000,1);
            frSTD_LC(i,cnt)=std(tmp);
            frMean_LC(i,cnt)=mean(tmp);
            frspont(i,cnt)=mean(FiringRate(200:1950,i,cnt));
            frSTDspont(i,cnt)=std(FiringRate(200:1950,i,cnt));
        end
    end
end
figure
subplot(2,1,1)
hold on
plot(kmean(1:2:end,onOff>0))
line([0 8],[0,0],'color','k')

subplot(2,1,2)
hold on
plot(kmean(1:2:end,onOff<0))
line([0 8],[0,0],'color','k')

kmeanON=kmean(1:2:end,onOff>0);
kmeanOFF=kmean(1:2:end,onOff<0);
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateAccum','kmeanON','kmeanOFF')
formula=1:2:2*8;
save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frspont', 'frSTDspont','frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','formula')
% 
% path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
% load([path2save,'quick'],'white_flash','black_flash','names')
% figure
% cc=1;
% for i=1:2:16
%     subplot(2,4,cc)
%     plot(mean(white_flash(:,i:i+1,26),2))
%     cc=cc+1;
% end

    

j=10
m=9
plot(correctedProtocols(:,1,m),correctedProtocols(:,2,m)+30,'k')
hold on
tmp=mean(FiringRate(:,m:m+1,j),2);
plot(tmp,'r')
legend('Stimulus','Firing Rate, ON cell')
for i=1:600:3600
    b=mean(tmp(correctedProtocols(i,1,m)+2000:correctedProtocols(i+600,1,m)-500))
    line([correctedProtocols(i,1,m)+2000 correctedProtocols(i+600,1,m)-500],[b,b],'color','b','linewidth',3)
    c=std(tmp(correctedProtocols(i,1,m)+2000:correctedProtocols(i+600,1,m)-500))
    line([correctedProtocols(i+600,1,m)-1000 correctedProtocols(i+600,1,m)-1000],[b-c,b+c],'color','b','linewidth',3)
end
for i=1:600:3601
    line([correctedProtocols(i,1,m) correctedProtocols(i,1,m)],[0 100],'color','g')
end
xlabel('Time')
ylabel('Firing Rate and Contrast')
title('Example of an ON cell with "wrong" behavior, ND4')

%% 20130301_1
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130301_1'
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

FiringRate=zeros(65000,size(SpikeCount,1),size(SpikeCount,2));
for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,65000);
            FiringRate(:,i,cnt)=tmp(121:end-120);
            tmp=FiringRate(:,i,cnt);
            for mc=1:6
                take=correctedProtocols(600*(mc-1)+1,1,i)+2000:correctedProtocols(600*mc,1,i)-500;
                tmp1(:,mc)=tmp(take(1:7000));
            end
            tmp=reshape(tmp1(:,1:2:end),21000,1);
            frSTD_HC(i,cnt)=std(tmp);
            frMean_HC(i,cnt)=mean(tmp);
            
            tmp=reshape(tmp1(:,2:2:end),21000,1);
            frSTD_LC(i,cnt)=std(tmp);
            frMean_LC(i,cnt)=mean(tmp);
            frspont(i,cnt)=mean(FiringRate(200:1950,i,cnt));
            frSTDspont(i,cnt)=std(FiringRate(200:1950,i,cnt));
        end
    end
end
figure
subplot(2,1,1)
hold on
plot(kmean(1:2:end,onOff>0))
line([0 8],[0,0],'color','k')

subplot(2,1,2)
hold on
plot(kmean(1:2:end,onOff<0))
line([0 8],[0,0],'color','k')


load('/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateAccum','kmeanON','kmeanOFF')
kmeanON=[kmeanON kmean(1:2:end,onOff>0)];
kmeanOFF=[kmeanOFF kmean(1:2:end,onOff<0)];
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateAccum','kmeanON','kmeanOFF')
formula=1:2:2*8;
save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frspont', 'frSTDspont','frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','formula')

%% 20130301_2
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130301_2'
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

FiringRate=zeros(65000,size(SpikeCount,1),size(SpikeCount,2));
for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,65000);
            FiringRate(:,i,cnt)=tmp(121:end-120);
            tmp=FiringRate(:,i,cnt);
            for mc=1:6
                take=correctedProtocols(600*(mc-1)+1,1,i)+2000:correctedProtocols(600*mc,1,i)-500;
                tmp1(:,mc)=tmp(take(1:7000));
            end
            tmp=reshape(tmp1(:,1:2:end),21000,1);
            frSTD_HC(i,cnt)=std(tmp);
            frMean_HC(i,cnt)=mean(tmp);
            
            tmp=reshape(tmp1(:,2:2:end),21000,1);
            frSTD_LC(i,cnt)=std(tmp);
            frMean_LC(i,cnt)=mean(tmp);
            frspont(i,cnt)=mean(FiringRate(200:1950,i,cnt));
            frSTDspont(i,cnt)=std(FiringRate(200:1950,i,cnt));
        end
    end
end
figure
subplot(2,1,1)
hold on
plot(kmean(1:2:end,onOff>0))
line([0 8],[0,0],'color','k')

subplot(2,1,2)
hold on
plot(kmean(1:2:end,onOff<0))
line([0 8],[0,0],'color','k')


load('/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateAccum','kmeanON','kmeanOFF')
kmeanON=[kmeanON kmean(1:2:end,onOff>0)];
kmeanOFF=[kmeanOFF kmean(1:2:end,onOff<0)];
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateAccum','kmeanON','kmeanOFF')
formula=1:2:2*8;
save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frspont', 'frSTDspont','frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','formula')

%% 20121023
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20121023'
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

FiringRate=zeros(65000,size(SpikeCount,1),size(SpikeCount,2));
for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,65000);
            FiringRate(:,i,cnt)=tmp(121:end-120);
        end
    end
end

figure
cnt=20
tmp=FiringRate(:,97:6:96+24,cnt);
plot(tmp)



for cnt=1:size(FiringRate,3)
i=1;
    for km=1:24:24*8
        tmp=FiringRate(:,km:2:(km+23),cnt);
        tmp=reshape(tmp(2501:62500,:),720000,1);
        
        frSTD_HC(i,cnt)=std(tmp);
        frMean_HC(i,cnt)=mean(tmp);
        
        tmp=FiringRate(:,(km+1):2:(km+23),cnt);
        tmp=reshape(tmp(2501:62500,:),720000,1);
        frSTD_LC(i,cnt)=std(tmp);
        frMean_LC(i,cnt)=mean(tmp);
        a=FiringRate(200:1950,km:km+23,cnt);
        frspont(i,cnt)=mean(a(:));
        frSTDspont(i,cnt)=std(a(:));
        i=i+1;
    end
end

figure
subplot(2,1,1)
plot(kmean(:,onOff>0))
line([0 8],[0,0],'color','k')

subplot(2,1,2)
plot(kmean(:,onOff<0))
line([0 8],[0,0],'color','k')


load('/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateAccum','kmeanON','kmeanOFF')
kmeanON=[kmeanON kmean(:,onOff>0)];
kmeanOFF=[kmeanOFF kmean(:,onOff<0)];
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateAccum','kmeanON','kmeanOFF')
formula=1:8;
save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frspont', 'frSTDspont','frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','formula')
% figure
% subplot(2,1,1)
% plot(kstd(:,onOff>0))
% line([0 8],[0,0],'color','k')
% 
% subplot(2,1,2)
% plot(kstd(:,onOff<0))
% line([0 8],[0,0],'color','k')


% 
% path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
% load([path2save,'quick'],'white_flash','black_flash','names')
% figure
% cc=1;
% for i=1:8
%     subplot(2,4,cc)
%     plot(white_flash(:,i,6))
%     cc=cc+1;
% end


%% 20130302
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130302'
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

FiringRate=zeros(65000,size(SpikeCount,1),size(SpikeCount,2));
for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,65000);
            FiringRate(:,i,cnt)=tmp(121:end-120);
            tmp=FiringRate(:,i,cnt);
            for mc=1:6
                take=correctedProtocols(600*(mc-1)+1,1,i)+2000:correctedProtocols(600*mc,1,i)-500;
                tmp1(:,mc)=tmp(take(1:7000));
            end
            tmp=reshape(tmp1(:,1:2:end),21000,1);
            frSTD_HC(i,cnt)=std(tmp);
            frMean_HC(i,cnt)=mean(tmp);
            
            tmp=reshape(tmp1(:,2:2:end),21000,1);
            frSTD_LC(i,cnt)=std(tmp);
            frMean_LC(i,cnt)=mean(tmp);
            frspont(i,cnt)=mean(FiringRate(200:1950,i,cnt));
            frSTDspont(i,cnt)=std(FiringRate(200:1950,i,cnt));
        end
    end
end
figure
subplot(2,1,1)
hold on
plot(kmean(1:3:end,onOff>0))
line([0 8],[0,0],'color','k')

subplot(2,1,2)
hold on
plot(kmean(1:3:end,onOff<0))
line([0 8],[0,0],'color','k')

load('/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateAccum','kmeanON','kmeanOFF')
kmeanON=[kmeanON kmean(1:3:end,onOff>0)];
kmeanOFF=[kmeanOFF kmean(1:3:end,onOff<0)];
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateAccum','kmeanON','kmeanOFF')
formula=1:3:3*8;
save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frspont', 'frSTDspont','frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','formula')


%% 20130302_1
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130302_1'
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

FiringRate=zeros(65000,size(SpikeCount,1),size(SpikeCount,2));
for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,65000);
            FiringRate(:,i,cnt)=tmp(121:end-120);
            tmp=FiringRate(:,i,cnt);
            for mc=1:6
                take=correctedProtocols(600*(mc-1)+1,1,i)+2000:correctedProtocols(600*mc,1,i)-500;
                tmp1(:,mc)=tmp(take(1:7000));
            end
            tmp=reshape(tmp1(:,1:2:end),21000,1);
            frSTD_HC(i,cnt)=std(tmp);
            frMean_HC(i,cnt)=mean(tmp);
            
            tmp=reshape(tmp1(:,2:2:end),21000,1);
            frSTD_LC(i,cnt)=std(tmp);
            frMean_LC(i,cnt)=mean(tmp);
            frspont(i,cnt)=mean(FiringRate(200:1950,i,cnt));
            frSTDspont(i,cnt)=std(FiringRate(200:1950,i,cnt));
        end
    end
end
figure
subplot(2,1,1)
hold on
plot(kmean(1:3:end,onOff>0))
line([0 8],[0,0],'color','k')

subplot(2,1,2)
hold on
plot(kmean(1:3:end,onOff<0))
line([0 8],[0,0],'color','k')

load('/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateAccum','kmeanON','kmeanOFF')
kmeanON=[kmeanON kmean(1:3:end,onOff>0)];
kmeanOFF=[kmeanOFF kmean(1:3:end,onOff<0)];
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateAccum','kmeanON','kmeanOFF')
formula=1:3:3*8;
save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frspont', 'frSTDspont','frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','formula')

%% 20130227
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20130227'
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

FiringRate=zeros(65000,size(SpikeCount,1),size(SpikeCount,2));
for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,65000);
            FiringRate(:,i,cnt)=tmp(121:end-120);
            tmp=FiringRate(:,i,cnt);
            for mc=1:6
                take=correctedProtocols(600*(mc-1)+1,1,i)+2000:correctedProtocols(600*mc,1,i)-500;
                tmp1(:,mc)=tmp(take(1:7000));
            end
            tmp=reshape(tmp1(:,1:2:end),21000,1);
            frSTD_HC(i,cnt)=std(tmp);
            frMean_HC(i,cnt)=mean(tmp);
            
            tmp=reshape(tmp1(:,2:2:end),21000,1);
            frSTD_LC(i,cnt)=std(tmp);
            frMean_LC(i,cnt)=mean(tmp);
            frspont(i,cnt)=mean(FiringRate(200:1950,i,cnt));
            frSTDspont(i,cnt)=std(FiringRate(200:1950,i,cnt));
        end
    end
end
figure
subplot(2,1,1)
hold on
plot(kmean(5:4:end,onOff>0))
line([0 8],[0,0],'color','k')

subplot(2,1,2)
hold on
plot(kmean(5:4:end,onOff<0))
line([0 8],[0,0],'color','k')

load('/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateAccum','kmeanON','kmeanOFF')
kmeanON=[kmeanON kmean(5:4:end,onOff>0)];
kmeanOFF=[kmeanOFF kmean(5:4:end,onOff<0)];
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateAccum','kmeanON','kmeanOFF')
formula=5:4:8*4+1;
save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frspont', 'frSTDspont','frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','formula')


%% 20120329
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20120329'
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

FiringRate=zeros(310000,size(SpikeCount,1),size(SpikeCount,2));
for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,310000);
            FiringRate(:,i,cnt)=tmp(121:end-120);
            tmp=FiringRate(:,i,cnt);
            for mc=1:10
                take=correctedProtocols(1800*(mc-1)+1,1,i)+2000:correctedProtocols(1800*mc,1,i)-500;
                tmp1(:,mc)=tmp(take(1:27500));
            end
            tmp=reshape(tmp1(:,1:2:end),27500*5,1);
            frSTD_HC(i,cnt)=std(tmp);
            frMean_HC(i,cnt)=mean(tmp);
            
            tmp=reshape(tmp1(:,2:2:end),27500*5,1);
            frSTD_LC(i,cnt)=std(tmp);
            frMean_LC(i,cnt)=mean(tmp);
            frspont(i,cnt)=mean(FiringRate(200:1950,i,cnt));
            frSTDspont(i,cnt)=std(FiringRate(200:1950,i,cnt));
        end
    end
end
figure
subplot(2,1,1)
hold on
plot(kmean(1:4:32,onOff>0))
line([0 8],[0,0],'color','k')

subplot(2,1,2)
hold on
plot(kmean(1:4:32,onOff<0))
line([0 8],[0,0],'color','k')


load('/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateAccum','kmeanON','kmeanOFF')
kmeanON=[kmeanON kmean(1:4:32,onOff>0)];
kmeanOFF=[kmeanOFF kmean(1:4:32,onOff<0)];
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateAccum','kmeanON','kmeanOFF')

formula=1:4:32;
save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frspont', 'frSTDspont','frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','formula')

%% 20120627
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20120627'
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

FiringRate=zeros(250000,size(SpikeCount,1),size(SpikeCount,2));
for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,250000);
            FiringRate(:,i,cnt)=tmp(121:end-120);
            tmp=FiringRate(:,i,cnt);
            for mc=1:8
                take=correctedProtocols(1800*(mc-1)+1,1,i)+2000:correctedProtocols(1800*mc,1,i)-500;
                tmp1(:,mc)=tmp(take(1:27500));
            end
            
            tmp=reshape(tmp1(:,1:2:end),27500*4,1);
            frSTD_HC(i,cnt)=std(tmp);
            frMean_HC(i,cnt)=mean(tmp);
            
            tmp=reshape(tmp1(:,2:2:end),27500*4,1);
            frSTD_LC(i,cnt)=std(tmp);
            frMean_LC(i,cnt)=mean(tmp);
            frspont(i,cnt)=mean(FiringRate(200:1950,i,cnt));
            frSTDspont(i,cnt)=std(FiringRate(200:1950,i,cnt));
        end
    end
end
figure
subplot(2,1,1)
hold on
plot(kmean(1:6:48,onOff>0))
line([0 8],[0,0],'color','k')

subplot(2,1,2)
hold on
plot(kmean(1:6:48,onOff<0))
line([0 8],[0,0],'color','k')


load('/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateAccum','kmeanON','kmeanOFF')
kmeanON=[kmeanON kmean(1:6:48,onOff>0)];
kmeanOFF=[kmeanOFF kmean(1:6:48,onOff<0)];
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateAccum','kmeanON','kmeanOFF')
formula=1:6:48;
save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frspont', 'frSTDspont','frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','formula')

%% 20120714
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20120714'
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

FiringRate=zeros(250000,size(SpikeCount,1)+6,size(SpikeCount,2));
for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,250000);
            FiringRate(:,i+6,cnt)=tmp(121:end-120);
            tmp=FiringRate(:,i+6,cnt);
            for mc=1:8
                take=correctedProtocols(1800*(mc-1)+1,1,i)+2000:correctedProtocols(1800*mc,1,i)-500;
                tmp1(:,mc)=tmp(take(1:27500));
            end
            
            tmp=reshape(tmp1(:,1:2:end),27500*4,1);
            frSTD_HC(i+6,cnt)=std(tmp);
            frMean_HC(i+6,cnt)=mean(tmp);
            
            tmp=reshape(tmp1(:,2:2:end),27500*4,1);
            frSTD_LC(i+6,cnt)=std(tmp);
            frMean_LC(i+6,cnt)=mean(tmp);
            frspont(i+6,cnt)=mean(FiringRate(200:1950,i+6,cnt));
            frSTDspont(i+6,cnt)=std(FiringRate(200:1950,i+6,cnt));
        end
    end
end
formula=1:6:48;
save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frspont', 'frSTDspont','frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','formula')


%% 20120902_1 and 20120902_2
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20120902_2'
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

FiringRate=zeros(310000,size(SpikeCount,1),size(SpikeCount,2));
for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,310000);
            FiringRate(:,i,cnt)=tmp(121:end-120);
            tmp=FiringRate(:,i,cnt);
            for mc=1:10
                take=correctedProtocols(1800*(mc-1)+1,1,i)+2000:correctedProtocols(1800*mc,1,i)-500;
                tmp1(:,mc)=tmp(take(1:27500));
            end
            tmp=reshape(tmp1(:,1:2:end),27500*5,1);
            frSTD_HC(i,cnt)=std(tmp);
            frMean_HC(i,cnt)=mean(tmp);
            
            tmp=reshape(tmp1(:,2:2:end),27500*5,1);
            frSTD_LC(i,cnt)=std(tmp);
            frMean_LC(i,cnt)=mean(tmp);
            frspont(i,cnt)=mean(FiringRate(200:1950,i,cnt));
            frSTDspont(i,cnt)=std(FiringRate(200:1950,i,cnt));
        end
    end
end


formula=1:4:12;
save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frspont', 'frSTDspont','frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','formula')

%% 20121023_1 and 20121026_1
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20121026_1'
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])


FiringRate=zeros(65000,size(SpikeCount,1),size(SpikeCount,2));
for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,65000);
            FiringRate(:,i,cnt)=tmp(121:end-120);
        end
    end
end


for cnt=1:size(FiringRate,3)
i=1;
    for km=1:24:24*5
        tmp=FiringRate(:,km:2:(km+23),cnt);
        tmp=reshape(tmp(2501:62500,:),720000,1);
        
        frSTD_HC(i,cnt)=std(tmp);
        frMean_HC(i,cnt)=mean(tmp);
        
        tmp=FiringRate(:,(km+1):2:(km+23),cnt);
        tmp=reshape(tmp(2501:62500,:),720000,1);
        frSTD_LC(i,cnt)=std(tmp);
        frMean_LC(i,cnt)=mean(tmp);
        a=FiringRate(200:1950,km:km+23,cnt);
        frspont(i,cnt)=mean(a(:));
        frSTDspont(i,cnt)=std(a(:));
        i=i+1;
    end
end


formula=1:5;
save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frspont', 'frSTDspont','frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','formula')
