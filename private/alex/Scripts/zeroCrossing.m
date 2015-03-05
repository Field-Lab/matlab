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
        for k=0:3
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
        for k=0:3
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
save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')

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
        for k=0:3
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
        for k=0:3
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
save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')

%% 20121023
clear
boundaries=[155 130 115 85 75 75 85 85];
date='20121023';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

cnt=1;
for i=1:24:24*8
    for j=1:size(LinearFilter,4)
        t=reshape(LinearFilter(:,1,i:2:i+23,j),500,12);
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
        
        t=reshape(LinearFilter(:,1,i+1:2:i+23,j),500,12);
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

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')

%% 20121023_1
clear
boundaries=[155 130 115 85 75 75 85 85];
date='20121023_1';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

cnt=1;
for i=1:24:24*5
    for j=1:size(LinearFilter,4)
        t=reshape(LinearFilter(:,1,i:2:i+23,j),500,12);
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
        
        t=reshape(LinearFilter(:,1,i+1:2:i+23,j),500,12);
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

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')

%% 20121026_1
clear
boundaries=[155 130 115 85 75 75 85 85];
date='20121026_1';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

cnt=1;
for i=1:24:24*5
    for j=1:size(LinearFilter,4)
        t=reshape(LinearFilter(:,1,i:2:i+23,j),500,12);
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
        
        t=reshape(LinearFilter(:,1,i+1:2:i+23,j),500,12);
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

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')


%% 20130220
clear
boundaries=[155 130 115 85 75 75 85 85];
date='20130220';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

cnt=1;
for i=11:9:82
    for j=1:size(LinearFilter,4)
        t=reshape(LinearFilter(:,1,i:i+8,j),500,9);
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

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')

%% 20130220_1 20130224 20130225 20130226 - just change the date 
clear
boundaries=[155 130 115 85 75 75 85 85];
date='20130226';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])
cnt=1;
for i=10:9:81
    for j=1:size(LinearFilter,4)
        t=reshape(LinearFilter(:,1,i:i+8,j),500,9);
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

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')


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
        for k=0:3
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
        for k=0:3
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
save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')


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
        for k=0:1
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
        for k=0:1
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

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')



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
        for k=0:1
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
        for k=0:1
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

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')




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
        for k=0:1
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
        for k=0:1
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

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')


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
        for k=0:2
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
        for k=0:2
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

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')


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
        for k=0:2
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
        for k=0:2
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

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')



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
        for k=0:5
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
        for k=0:5
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

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')


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
        for k=0:5
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
        for k=0:5
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

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')



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
        for k=0:2
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
        for k=0:2
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

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')



%% Plots

clear

dates=cell(7,1);
dates{1}='20130301'
dates{2}='20130301_1'
dates{3}='20130301_2'
dates{4}='20130302'
dates{5}='20130302_1'
dates{6}='20120329'
dates{7}='20121023'

zc_HC_ON=[];zc_LC_ON=[];
zc_HC_OFF=[];zc_LC_OFF=[];

for datesCNT=1:7
    date=dates{datesCNT}
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date],'onOff','zc_HC','zc_LC')
    zc_HC_ON=[zc_HC_ON; zc_HC(onOff>0,:)];
    zc_LC_ON=[zc_LC_ON; zc_LC(onOff>0,:)];
    zc_HC_OFF=[zc_HC_OFF; zc_HC(onOff<0,:)];
    zc_LC_OFF=[zc_LC_OFF; zc_LC(onOff<0,:)];
end




subplot(2,1,1)
plot((zc_HC_ON-zc_LC_ON)')
line([0 9],[0,0],'color','k','linewidth',2)
subplot(2,1,2)
plot((zc_HC_OFF-zc_LC_OFF)')
line([0 9],[0,0],'color','k','linewidth',2)

figure
for i=1:8
    subplot(4,2,i)
    plot(sort((zc_HC_ON(:,i)-zc_LC_ON(:,i))),'.')
    line([0 size(zc_HC_ON,1)],[0,0],'color','k','linewidth',2)
    axis([0 size(zc_HC_ON,1) -15 15])
end

figure
for i=1:8
    subplot(4,2,i)
    plot(sort((zc_HC_OFF(:,i)-zc_LC_OFF(:,i))),'.')
    line([0 size(zc_HC_OFF,1)],[0,0],'color','k','linewidth',2)
    axis([0 size(zc_HC_OFF,1) -15 15])
end


a=zc_HC_ON-zc_LC_ON;
k1=sum(a>=2&a<=10)/size(a,1)*100;
a=zc_HC_OFF-zc_LC_OFF;
k=sum(a>=2&a<=10)/size(a,1)*100;
figure
bar([k; k1]')
legend('OFF','ON')
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND')
ylabel('% of total')
title('portion of cells speeding up at low contrast')



a=zc_HC_ON-zc_LC_ON;
k1=sum(a<2&a>-2)/size(a,1)*100;
a=zc_HC_OFF-zc_LC_OFF;
k=sum(a<2&a>-2)/size(a,1)*100;
figure
bar([k; k1]')
legend('OFF','ON')
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND')
ylabel('% of total')
title('portion of cells not changing at low contrast')



a=zc_HC_ON-zc_LC_ON;
k1=sum(a<=-2&a>-10)/size(a,1)*100;
a=zc_HC_OFF-zc_LC_OFF;
k=sum(a<=-2&a>-10)/size(a,1)*100;
figure
bar([k; k1]')
legend('OFF','ON')
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND')
ylabel('% of total')
title('portion of cells slowing down at low contrast')










a=zc_HC_ON-zc_LC_ON;
k=sum(a>=-2&a<=2)/size(a,1)*100;
k1=sum(a>2&a<=10)/size(a,1)*100;
k2=sum(a<-2&a>-10)/size(a,1)*100;
figure
subplot(2,1,1)
bar([k1; k; k2]')
axis([0 9 0 50])
legend({'speed up','no change', 'slow down'})
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND')
ylabel('% of total')
title('ON cells: changes at low contrast')
a=zc_HC_OFF-zc_LC_OFF;
k=sum(a>=-2&a<=2)/size(a,1)*100;
k1=sum(a>2&a<=10)/size(a,1)*100;
k2=sum(a<-2&a>-10)/size(a,1)*100;
subplot(2,1,2)
bar([k1; k; k2]')
axis([0 9 0 50])
legend({'speed up','no change', 'slow down'})
set(gca,'xtick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND')
ylabel('% of total')
title('OFF cells: changes at low contrast')





%% Example 
%% Plots for 20121023
date='20121023';
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

figure
set(gcf,'position',[560          17        1677         931])
nds='87654321';

for j=1:size(LinearFilter,4)
    cnt=1;
    for i=1:24:24*8
        subplot(2,4,cnt)
        hold off
        t=reshape(LinearFilter(:,1,i:2:i+23,j),500,12);
        t=t-repmat(mean(t),size(t,1),1);
        t=t./repmat(sum(abs(t)),size(t,1),1);
        t=mean(t,2);
        plot(t(1:350),'linewidth',2)
        hold on        
        t=reshape(LinearFilter(:,1,i+1:2:i+23,j),500,12);
        t=t-repmat(mean(t),size(t,1),1);
        t=t./repmat(sum(abs(t)),size(t,1),1);
        t=mean(t,2);
        plot(t(1:350),'r','linewidth',2)
        line([0 350],[0,0],'color','k')
        axis([0,350 -0.02 0.02])
        title(['ND',nds(cnt)])
        cnt=cnt+1;
    end
    subplot('Position',[0.5 0.96 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([names{j}, '   Linear Filters High (blue) and Low (red) Contrast Normalized to mean 0, area 1'],'FontSize',12,'FontWeight','bold','Interpreter','None')
end

