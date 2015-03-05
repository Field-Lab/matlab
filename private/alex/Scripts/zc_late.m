boundaries=[155 130 115 85 75 75 85 85];

%% 20120902_1
clear
boundaries=[155 130 115 85 75 75 85 85];
date='20120902_1';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

zc_HC=[];
zc_LC=[];
for j=1:size(LinearFilter,4)
    cnt=1;
    for i=1:4:4*4

        m=[];
        for k=2:3
            t=LinearFilter(:,1:2:end,i+k,j);
            t=t-repmat(mean(t),size(t,1),1);
            t=t./repmat(sum(abs(t)),size(t,1),1);
            m=[m t];
        end
        m=nanmean(m,2);
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_HC(j,cnt)=myFit(m(zc_HC(j,cnt)-boundaries(cnt)+1:zc_HC(j,cnt)),zc_HC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_HC(j,cnt)=0;
            zc_HC(j,cnt)=0;
        end
        m=[];
        for k=2:3
            t=LinearFilter(:,2:2:end,i+k,j);
            t=t-repmat(mean(t),size(t,1),1);
            t=t./repmat(sum(abs(t)),size(t,1),1);
            m=[m t];
        end
        m=nanmean(m,2);
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_LC(j,cnt)=myFit(m(zc_LC(j,cnt)-boundaries(cnt)+1:zc_LC(j,cnt)),zc_LC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_LC(j,cnt)=0;
            zc_LC(j,cnt)=0;
        end
        cnt=cnt+1;
    end
end
save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_1',date],'zc_HC','zc_LC','peak_HC','peak_LC')

%% 20120902_1
clear
boundaries=[155 130 115 85 75 75 85 85];
date='20120902_2';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

zc_HC=[];
zc_LC=[];
for j=1:size(LinearFilter,4)
    cnt=1;
    for i=1:4:4*3

        m=[];
        for k=2:3
            t=LinearFilter(:,1:2:end,i+k,j);
            t=t-repmat(mean(t),size(t,1),1);
            t=t./repmat(sum(abs(t)),size(t,1),1);
            m=[m t];
        end
        m=nanmean(m,2);
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_HC(j,cnt)=myFit(m(zc_HC(j,cnt)-boundaries(cnt)+1:zc_HC(j,cnt)),zc_HC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_HC(j,cnt)=0;
            zc_HC(j,cnt)=0;
        end
        m=[];
        for k=2:3
            t=LinearFilter(:,2:2:end,i+k,j);
            t=t-repmat(mean(t),size(t,1),1);
            t=t./repmat(sum(abs(t)),size(t,1),1);
            m=[m t];
        end
        m=nanmean(m,2);
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_LC(j,cnt)=myFit(m(zc_LC(j,cnt)-boundaries(cnt)+1:zc_LC(j,cnt)),zc_LC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_LC(j,cnt)=0;
            zc_LC(j,cnt)=0;
        end
        cnt=cnt+1;
    end
end
save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_1',date],'zc_HC','zc_LC','peak_HC','peak_LC')

%% 20121023
clear
boundaries=[155 130 115 85 75 75 85 85];
date='20121023';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

cnt=1;
for i=1:24:24*8
    for j=1:size(LinearFilter,4)
        t=reshape(LinearFilter(:,1,i+16:2:i+23,j),500,4);
        t=t-repmat(mean(t),size(t,1),1);
        t=t./repmat(sum(abs(t)),size(t,1),1);
        m=nanmean(t,2);
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_HC(j,cnt)=myFit(m(zc_HC(j,cnt)-boundaries(cnt)+1:zc_HC(j,cnt)),zc_HC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_HC(j,cnt)=0;
            zc_HC(j,cnt)=0;
        end
        
        t=reshape(LinearFilter(:,1,i+17:2:i+23,j),500,4);
        t=t-repmat(mean(t),size(t,1),1);
        t=t./repmat(sum(abs(t)),size(t,1),1);
        m=nanmean(t,2);
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_LC(j,cnt)=myFit(m(zc_LC(j,cnt)-boundaries(cnt)+1:zc_LC(j,cnt)),zc_LC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_LC(j,cnt)=0;
            zc_LC(j,cnt)=0;
        end
    end
    cnt=cnt+1;
end

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_1',date],'zc_HC','zc_LC','peak_HC','peak_LC')

%% 20121023_1
clear
boundaries=[155 130 115 85 75 75 85 85];
date='20121023_1';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

cnt=1;
for i=1:24:24*5
    for j=1:size(LinearFilter,4)
        t=reshape(LinearFilter(:,1,i+16:2:i+23,j),500,4);
        t=t-repmat(mean(t),size(t,1),1);
        t=t./repmat(sum(abs(t)),size(t,1),1);
        m=nanmean(t,2);
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_HC(j,cnt)=myFit(m(zc_HC(j,cnt)-boundaries(cnt)+1:zc_HC(j,cnt)),zc_HC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_HC(j,cnt)=0;
            zc_HC(j,cnt)=0;
        end
        
        t=reshape(LinearFilter(:,1,i+17:2:i+23,j),500,4);
        t=t-repmat(mean(t),size(t,1),1);
        t=t./repmat(sum(abs(t)),size(t,1),1);
        m=nanmean(t,2);
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_LC(j,cnt)=myFit(m(zc_LC(j,cnt)-boundaries(cnt)+1:zc_LC(j,cnt)),zc_LC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_LC(j,cnt)=0;
            zc_LC(j,cnt)=0;
        end
    end
    cnt=cnt+1;
end

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_1',date],'zc_HC','zc_LC','peak_HC','peak_LC')

%% 20121026_1
clear
boundaries=[155 130 115 85 75 75 85 85];
date='20121026_1';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

cnt=1;
for i=1:24:24*5
    for j=1:size(LinearFilter,4)
        t=reshape(LinearFilter(:,1,i+16:2:i+23,j),500,4);
        t=t-repmat(mean(t),size(t,1),1);
        t=t./repmat(sum(abs(t)),size(t,1),1);
        m=nanmean(t,2);
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_HC(j,cnt)=myFit(m(zc_HC(j,cnt)-boundaries(cnt)+1:zc_HC(j,cnt)),zc_HC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_HC(j,cnt)=0;
            zc_HC(j,cnt)=0;
        end
        
        t=reshape(LinearFilter(:,1,i+17:2:i+23,j),500,4);
        t=t-repmat(mean(t),size(t,1),1);
        t=t./repmat(sum(abs(t)),size(t,1),1);
        m=nanmean(t,2);
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_LC(j,cnt)=myFit(m(zc_LC(j,cnt)-boundaries(cnt)+1:zc_LC(j,cnt)),zc_LC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_LC(j,cnt)=0;
            zc_LC(j,cnt)=0;
        end
    end
    cnt=cnt+1;
end

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_1',date],'zc_HC','zc_LC','peak_HC','peak_LC')


%% 20130220
clear
boundaries=[155 130 115 85 75 75 85 85];
date='20130220';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

cnt=1;
for i=11:9:82
    for j=1:size(LinearFilter,4)
        t=reshape(LinearFilter(:,1,i+5:i+8,j),500,4);
        t=t-repmat(mean(t),size(t,1),1);
        t=t./repmat(sum(abs(t)),size(t,1),1);
        m=nanmean(t,2);
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_HC(j,cnt)=myFit(m(zc_HC(j,cnt)-boundaries(cnt)+1:zc_HC(j,cnt)),zc_HC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_HC(j,cnt)=0;
            zc_HC(j,cnt)=0;
        end     
        zc_LC(j,cnt)=0;
        peak_LC(j,cnt)=0;
    end
    cnt=cnt+1;
end

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_1',date],'zc_HC','zc_LC','peak_HC','peak_LC')

%% 20130220_1 20130224 20130225 20130226 - just change the date 
clear
boundaries=[155 130 115 85 75 75 85 85];
date='20130226';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])
cnt=1;
for i=10:9:81
    for j=1:size(LinearFilter,4)
        t=reshape(LinearFilter(:,1,i+5:i+8,j),500,4);
        t=t-repmat(mean(t),size(t,1),1);
        t=t./repmat(sum(abs(t)),size(t,1),1);
        m=nanmean(t,2);
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_HC(j,cnt)=myFit(m(zc_HC(j,cnt)-boundaries(cnt)+1:zc_HC(j,cnt)),zc_HC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_HC(j,cnt)=0;
            zc_HC(j,cnt)=0;
        end
        zc_LC(j,cnt)=0;
        peak_LC(j,cnt)=0;
    end
    cnt=cnt+1;
end

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_1',date],'zc_HC','zc_LC','peak_HC','peak_LC')


%% 20120329
clear
boundaries=[155 130 115 85 75 75 85 85];
date='20120329';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

zc_HC=[];
zc_LC=[];
for j=1:size(LinearFilter,4)
    cnt=1;
    for i=1:4:8*4

        m=[];
        for k=2:3
            t=LinearFilter(:,1:2:end,i+k,j);
            t=t-repmat(mean(t),size(t,1),1);
            t=t./repmat(sum(abs(t)),size(t,1),1);
            m=[m t];
        end
        m=nanmean(m,2);
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_HC(j,cnt)=myFit(m(zc_HC(j,cnt)-boundaries(cnt)+1:zc_HC(j,cnt)),zc_HC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_HC(j,cnt)=0;
            zc_HC(j,cnt)=0;
        end
        m=[];
        for k=2:3
            t=LinearFilter(:,2:2:end,i+k,j);
            t=t-repmat(mean(t),size(t,1),1);
            t=t./repmat(sum(abs(t)),size(t,1),1);
            m=[m t];
        end
        m=nanmean(m,2);
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_LC(j,cnt)=myFit(m(zc_LC(j,cnt)-boundaries(cnt)+1:zc_LC(j,cnt)),zc_LC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_LC(j,cnt)=0;
            zc_LC(j,cnt)=0;
        end
        cnt=cnt+1;
    end
end
save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_1',date],'zc_HC','zc_LC','peak_HC','peak_LC')


%% 20130301
clear
boundaries=[155 130 115 85 75 75 85 85];
date='20130301';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

zc_HC=[];
zc_LC=[];
for j=1:size(LinearFilter,4)
    cnt=1;
    for i=1:2:8*2

        m=[];
        for k=1
            t=LinearFilter(:,1:2:end,i+k,j);
            t=t-repmat(mean(t),size(t,1),1);
            t=t./repmat(sum(abs(t)),size(t,1),1);
            m=[m t];
        end
        m=nanmean(m,2);
        
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_HC(j,cnt)=myFit(m(zc_HC(j,cnt)-boundaries(cnt)+1:zc_HC(j,cnt)),zc_HC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_HC(j,cnt)=0;
            zc_HC(j,cnt)=0;
        end
        
        m=[];
        for k=1
            t=LinearFilter(:,2:2:end,i+k,j);
            t=t-repmat(mean(t),size(t,1),1);
            t=t./repmat(sum(abs(t)),size(t,1),1);
            m=[m t];
        end
        m=nanmean(m,2);
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_LC(j,cnt)=myFit(m(zc_LC(j,cnt)-boundaries(cnt)+1:zc_LC(j,cnt)),zc_LC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_LC(j,cnt)=0;
            zc_LC(j,cnt)=0;
        end
        cnt=cnt+1;
    end
end

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_1',date],'zc_HC','zc_LC','peak_HC','peak_LC')



%% 20130301_1
clear
boundaries=[155 130 115 85 75 75 85 85];
date='20130301_1';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

zc_HC=[];
zc_LC=[];
for j=1:size(LinearFilter,4)
    cnt=1;
    for i=1:2:8*2

        m=[];
        for k=1
            t=LinearFilter(:,1:2:end,i+k,j);
            t=t-repmat(mean(t),size(t,1),1);
            t=t./repmat(sum(abs(t)),size(t,1),1);
            m=[m t];
        end
        m=nanmean(m,2);
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_HC(j,cnt)=myFit(m(zc_HC(j,cnt)-boundaries(cnt)+1:zc_HC(j,cnt)),zc_HC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_HC(j,cnt)=0;
            zc_HC(j,cnt)=0;
        end
        
        
        m=[];
        for k=1
            t=LinearFilter(:,2:2:end,i+k,j);
            t=t-repmat(mean(t),size(t,1),1);
            t=t./repmat(sum(abs(t)),size(t,1),1);
            m=[m t];
        end
        m=nanmean(m,2);
        
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_LC(j,cnt)=myFit(m(zc_LC(j,cnt)-boundaries(cnt)+1:zc_LC(j,cnt)),zc_LC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_LC(j,cnt)=0;
            zc_LC(j,cnt)=0;
        end
        
        cnt=cnt+1;
    end
end

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_1',date],'zc_HC','zc_LC','peak_HC','peak_LC')




%% 20130301_2
clear
boundaries=[155 130 115 85 75 75 85 85];
date='20130301_2';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

zc_HC=[];
zc_LC=[];
for j=1:size(LinearFilter,4)
    cnt=1;
    for i=1:2:8*2

        m=[];
        for k=1
            t=LinearFilter(:,1:2:end,i+k,j);
            t=t-repmat(mean(t),size(t,1),1);
            t=t./repmat(sum(abs(t)),size(t,1),1);
            m=[m t];
        end
        m=nanmean(m,2);
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_HC(j,cnt)=myFit(m(zc_HC(j,cnt)-boundaries(cnt)+1:zc_HC(j,cnt)),zc_HC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_HC(j,cnt)=0;
            zc_HC(j,cnt)=0;
        end
        
        
        m=[];
        for k=1
            t=LinearFilter(:,2:2:end,i+k,j);
            t=t-repmat(mean(t),size(t,1),1);
            t=t./repmat(sum(abs(t)),size(t,1),1);
            m=[m t];
        end
        m=nanmean(m,2);
        
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_LC(j,cnt)=myFit(m(zc_LC(j,cnt)-boundaries(cnt)+1:zc_LC(j,cnt)),zc_LC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_LC(j,cnt)=0;
            zc_LC(j,cnt)=0;
        end
        cnt=cnt+1;
    end
end

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_1',date],'zc_HC','zc_LC','peak_HC','peak_LC')


%% 20130302
clear
boundaries=[155 130 115 85 75 75 85 85];
date='20130302';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

zc_HC=[];
zc_LC=[];
for j=1:size(LinearFilter,4)
    cnt=1;
    for i=1:3:8*3

        m=[];
        for k=2
            t=LinearFilter(:,1:2:end,i+k,j);
            t=t-repmat(mean(t),size(t,1),1);
            t=t./repmat(sum(abs(t)),size(t,1),1);
            m=[m t];
        end
        m=nanmean(m,2);
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_HC(j,cnt)=myFit(m(zc_HC(j,cnt)-boundaries(cnt)+1:zc_HC(j,cnt)),zc_HC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_HC(j,cnt)=0;
            zc_HC(j,cnt)=0;
        end
        
        m=[];
        for k=2
            t=LinearFilter(:,2:2:end,i+k,j);
            t=t-repmat(mean(t),size(t,1),1);
            t=t./repmat(sum(abs(t)),size(t,1),1);
            m=[m t];
        end
        m=nanmean(m,2);
        
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_LC(j,cnt)=myFit(m(zc_LC(j,cnt)-boundaries(cnt)+1:zc_LC(j,cnt)),zc_LC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_LC(j,cnt)=0;
            zc_LC(j,cnt)=0;
        end
        cnt=cnt+1;
    end
end

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_1',date],'zc_HC','zc_LC','peak_HC','peak_LC')


%% 20130302_1
clear
boundaries=[155 130 115 85 75 75 85 85];
date='20130302_1';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

zc_HC=[];
zc_LC=[];
for j=1:size(LinearFilter,4)
    cnt=1;
    for i=1:3:8*3

        m=[];
        for k=2
            t=LinearFilter(:,1:2:end,i+k,j);
            t=t-repmat(mean(t),size(t,1),1);
            t=t./repmat(sum(abs(t)),size(t,1),1);
            m=[m t];
        end
        m=nanmean(m,2);
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_HC(j,cnt)=myFit(m(zc_HC(j,cnt)-boundaries(cnt)+1:zc_HC(j,cnt)),zc_HC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_HC(j,cnt)=0;
            zc_HC(j,cnt)=0;
        end
        
        m=[];
        for k=2
            t=LinearFilter(:,2:2:end,i+k,j);
            t=t-repmat(mean(t),size(t,1),1);
            t=t./repmat(sum(abs(t)),size(t,1),1);
            m=[m t];
        end
        m=nanmean(m,2);
        
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_LC(j,cnt)=myFit(m(zc_LC(j,cnt)-boundaries(cnt)+1:zc_LC(j,cnt)),zc_LC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_LC(j,cnt)=0;
            zc_LC(j,cnt)=0;
        end
        cnt=cnt+1;
    end
end

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_1',date],'zc_HC','zc_LC','peak_HC','peak_LC')



%% 20120627
clear
boundaries=[155 130 115 85 75 75 85 85];
date='20120627';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

zc_HC=[];
zc_LC=[];
for j=1:size(LinearFilter,4)
    cnt=1;
    for i=1:6:8*6

        m=[];
        for k=4:5
            t=LinearFilter(:,1:2:end,i+k,j);
            t=t-repmat(mean(t),size(t,1),1);
            t=t./repmat(sum(abs(t)),size(t,1),1);
            m=[m t];
        end
        m=nanmean(m,2);
        
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_HC(j,cnt)=myFit(m(zc_HC(j,cnt)-boundaries(cnt)+1:zc_HC(j,cnt)),zc_HC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_HC(j,cnt)=0;
            zc_HC(j,cnt)=0;
        end
        
        m=[];
        for k=4:5
            t=LinearFilter(:,2:2:end,i+k,j);
            t=t-repmat(mean(t),size(t,1),1);
            t=t./repmat(sum(abs(t)),size(t,1),1);
            m=[m t];
        end
        m=nanmean(m,2);
        
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_LC(j,cnt)=myFit(m(zc_LC(j,cnt)-boundaries(cnt)+1:zc_LC(j,cnt)),zc_LC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_LC(j,cnt)=0;
            zc_LC(j,cnt)=0;
        end
        cnt=cnt+1;
    end
end

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_1',date],'zc_HC','zc_LC','peak_HC','peak_LC')


%% 20120714
clear
boundaries=[155 130 115 85 75 75 85 85];
date='20120714';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

zc_HC=zeros(size(LinearFilter,4),8);
zc_LC=zeros(size(LinearFilter,4),8);
peak_HC=zeros(size(LinearFilter,4),8);
peak_LC=zeros(size(LinearFilter,4),8);
for j=1:size(LinearFilter,4)
    cnt=2;
    for i=1:6:7*6

        m=[];
        for k=4:5
            t=LinearFilter(:,1:2:end,i+k,j);
            t=t-repmat(mean(t),size(t,1),1);
            t=t./repmat(sum(abs(t)),size(t,1),1);
            m=[m t];
        end
        m=nanmean(m,2);
        
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_HC(j,cnt)=myFit(m(zc_HC(j,cnt)-boundaries(cnt)+1:zc_HC(j,cnt)),zc_HC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_HC(j,cnt)=0;
            zc_HC(j,cnt)=0;
        end
        
        m=[];
        for k=4:5
            t=LinearFilter(:,2:2:end,i+k,j);
            t=t-repmat(mean(t),size(t,1),1);
            t=t./repmat(sum(abs(t)),size(t,1),1);
            m=[m t];
        end
        m=nanmean(m,2);
        
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_LC(j,cnt)=myFit(m(zc_LC(j,cnt)-boundaries(cnt)+1:zc_LC(j,cnt)),zc_LC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_LC(j,cnt)=0;
            zc_LC(j,cnt)=0;
        end
        cnt=cnt+1;
    end
end

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_1',date],'zc_HC','zc_LC','peak_HC','peak_LC')



%% 20130227
clear
boundaries=[155 130 115 85 75 75 85 85];
date='20130227';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

zc_HC=[];
zc_LC=[];
for j=1:size(LinearFilter,4)
    cnt=1;
    for i=1:3:8*3

        m=[];
        for k=2
            t=LinearFilter(:,1:2:end,i+k,j);
            t=t-repmat(mean(t),size(t,1),1);
            t=t./repmat(sum(abs(t)),size(t,1),1);
            m=[m t];
        end
        m=nanmean(m,2);
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_HC(j,cnt)=myFit(m(zc_HC(j,cnt)-boundaries(cnt)+1:zc_HC(j,cnt)),zc_HC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_HC(j,cnt)=0;
            zc_HC(j,cnt)=0;
        end
        
        
        m=[];
        for k=2
            t=LinearFilter(:,2:2:end,i+k,j);
            t=t-repmat(mean(t),size(t,1),1);
            t=t./repmat(sum(abs(t)),size(t,1),1);
            m=[m t];
        end
        m=nanmean(m,2);
        
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_LC(j,cnt)=myFit(m(zc_LC(j,cnt)-boundaries(cnt)+1:zc_LC(j,cnt)),zc_LC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_LC(j,cnt)=0;
            zc_LC(j,cnt)=0;
        end
        cnt=cnt+1;
    end
end

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_1',date],'zc_HC','zc_LC','peak_HC','peak_LC')
