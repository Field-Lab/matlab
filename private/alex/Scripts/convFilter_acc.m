cd('/mnt/muench_data/user/alexandra/scripts')

clear

dates=cell(15,1);
dates{1}='20130220';
dates{2}='20130220_1';
dates{3}='20130224';
dates{4}='20130225';
dates{5}='20130226';
dates{6}='20130227';
dates{7}='20130301';
dates{8}='20130301_1';
dates{9}='20130301_2';
dates{10}='20130302';
dates{11}='20130302_1';
dates{12}='20120329';
dates{13}='20120627';
dates{14}='20120714';
dates{15}='20121023';

ratio_area_all=cell(1,8);
peak_neg_all=cell(1,8);
peak_pos_all=cell(1,8);
ratio_area_all_LC=cell(1,8);
peak_neg_all_LC=cell(1,8);
peak_pos_all_LC=cell(1,8);
basFR=cell(1,8);
HC_pos=zeros(700,8,520);
HC_neg=zeros(700,8,520);
LC_pos=zeros(700,8,520);
LC_neg=zeros(700,8,520);
last_tt=0;
for datesCNT=1:15
    
    date=dates{datesCNT}
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date],'LinearFilter','path2save','dim','names','onOff')
    load([path2save,'convFilters_std'])
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frspont','frMean_HC','frMean_LC')
    
    clear peak_pos peak_neg ratio_peaks ratio_area
    for cnt=1:size(posFilter_std,3)
        for i=1:size(posFilter_std,2)
            if onOff(cnt)>0
                [k,m]=max(posFilter_std(201:450,i,cnt));
                [p,n]=max(negFilter_std(201:450,i,cnt));
            elseif onOff(cnt)<0
                [k,m]=min(posFilter_std(201:450,i,cnt));
                [p,n]=min(negFilter_std(201:450,i,cnt));
            end
            a=sum(abs((posFilter_std(201:450,i,cnt))));
            b=sum(abs((negFilter_std(201:450,i,cnt))));
            
            if onOff(cnt)~=0
                peak_pos(i,cnt)=m;
                peak_neg(i,cnt)=n;
                ratio_peaks(i,cnt)=(abs(k)-abs(p))/(abs(k)+abs(p));
                ratio_area(i,cnt)=(abs(a)-abs(b))/(abs(a)+abs(b));
            else
                peak_pos(i,cnt)=-1;
                peak_neg(i,cnt)=-1;
                ratio_peaks(i,cnt)=-1;
                ratio_area(i,cnt)=-1;
            end
        end
    end
    for nd=1:8
        high_contr=metastat(1:size(posFilter_std,2),1)==nd&metastat(1:size(posFilter_std,2),2)==0;
        
        
        tmp=mean(posFilter_std(:,high_contr,:),2);
        ttt=size(tmp,3);
        HC_pos(:,9-nd,last_tt+1:last_tt+ttt)=tmp;
        tmp=mean(negFilter_std(:,high_contr,:),2);
        ttt=size(tmp,3);
        HC_neg(:,9-nd,last_tt+1:last_tt+ttt)=tmp;        
        
        k=nanmean(ratio_area(high_contr,:));
        ratio_area_all{9-nd}=[ratio_area_all{9-nd} k];
        k=nanmean(peak_pos(high_contr,:));
        peak_pos_all{9-nd}=[peak_pos_all{9-nd} k];
        k=nanmean(peak_neg(high_contr,:));
        peak_neg_all{9-nd}=[peak_neg_all{9-nd} k];
        k=nanmean(basalFR(metastat(metastat(1:size(posFilter_std,2),2)==0,1)==nd,:));
        basFR{9-nd}=[basFR{9-nd} k];
        
        low_contr=metastat(1:size(posFilter_std,2),1)==nd&metastat(1:size(posFilter_std,2),2)==1;
        if sum(low_contr)>0
            k=nanmean(ratio_area(low_contr,:));
            ratio_area_all_LC{9-nd}=[ratio_area_all_LC{9-nd} k];
            k=nanmean(peak_pos(low_contr,:));
            peak_pos_all_LC{9-nd}=[peak_pos_all_LC{9-nd} k];
            k=nanmean(peak_neg(low_contr,:));
            peak_neg_all_LC{9-nd}=[peak_neg_all_LC{9-nd} k];
            
            tmp=mean(posFilter_std(:,low_contr,:),2);
            ttt=size(tmp,3);
            LC_pos(:,9-nd,last_tt+1:last_tt+ttt)=tmp;
            tmp=mean(negFilter_std(:,low_contr,:),2);
            ttt=size(tmp,3);
            LC_neg(:,9-nd,last_tt+1:last_tt+ttt)=tmp;
            
        else
            k=nanmean(ratio_area(high_contr,:));
            ratio_area_all_LC{9-nd}=[ratio_area_all_LC{9-nd} k*0];
            peak_pos_all_LC{9-nd}=[peak_pos_all_LC{9-nd} k*0];
            peak_neg_all_LC{9-nd}=[peak_neg_all_LC{9-nd} k*0];            
            LC_pos(:,9-nd,last_tt+1:last_tt+ttt)=0;
            LC_neg(:,9-nd,last_tt+1:last_tt+ttt)=0;            
        end
    end
    last_tt=last_tt+ttt;
    
end

ratio_area_HC=zeros(520,8);
peak_neg_HC=zeros(520,8);
peak_pos_HC=zeros(520,8);
ratio_area_LC=zeros(520,8);
peak_neg_LC=zeros(520,8);
peak_pos_LC=zeros(520,8);
bas_FR=zeros(520,8);
for i=1:8
    ratio_area_HC(:,i)=ratio_area_all{i}';
    peak_neg_HC(:,i)=peak_neg_all{i}';
    peak_pos_HC(:,i)=peak_pos_all{i}';
    ratio_area_LC(:,i)=ratio_area_all_LC{i}';
    peak_neg_LC(:,i)=peak_neg_all_LC{i}';
    peak_pos_LC(:,i)=peak_pos_all_LC{i}';
    bas_FR(:,i)=basFR{i}';
end

load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all'])




figure
for cnt=481+(1:39)
    subplot(6,7,cnt-481)
    hold on
    cc=1;
    for i=3        
        a=(HC_pos(:,3,cnt));
        plot(a)
         a=(HC_neg(:,3,cnt));
        plot(a,'r')
    end
    line([0 700],[0 0],'color','k')
    axis tight
    title(int2str(cnt-481))
end




for i=1:8
    subplot(2,4,i)
    cond=onOff>0&bothContr&bas_FR(:,i)'>10;
    sum(cond)

end




figure
cond=onOff>0&bothContr;

a=ratio_area_HC(cond,3:5);

% b=ratio_area_LC(cond,i);

plot(a')
title('high contrast area ratio, ON cells')




figure
for i=1:8
    subplot(2,4,i)
    cond=onOff>0&bothContr&zc_LC(i,:)>50&zc_HC(i,:)>50&ratio_area_HC(:,i)'<0.8;
    
    a=ratio_area_HC(cond,i);
    
    b=ratio_area_LC(cond,i);

    plot(a,b,'.r')
%     hold on
% %    cond=onOff<0&bothContr&zc_LC(i,:)>50&zc_HC(i,:)>50&ratio_area_HC(:,i)'<0.3;
%     a=reshape(HC_neg(:,i,cond),700,sum(cond));
%     plot(a,'b')
%     title(int2str((sum(cond)/sum(predC))*100))
% xlabel('pos')
% ylabel('neg')
% % line([0,0],[-30 30])
% line([-30 30],[0,0])
% axis([-30 30 -30 30])
end






figure
for i=1:8
    subplot(2,4,i)
    cond=onOff>0&bothContr&zc_LC(i,:)>50&zc_HC(i,:)>50&ratio_area_HC(:,i)'<0.2;
    predC=onOff>0&bothContr&zc_LC(i,:)>50&zc_HC(i,:)>50;
    
%     names{cond}
    a=reshape(HC_pos(:,i,cond),700,sum(cond));
    plot(a,'r')
    hold on
%    cond=onOff<0&bothContr&zc_LC(i,:)>50&zc_HC(i,:)>50&ratio_area_HC(:,i)'<0.3;
    a=reshape(HC_neg(:,i,cond),700,sum(cond));
    plot(a,'b')
    title(int2str((sum(cond)/sum(predC))*100))
end


figure
for i=1:8
    subplot(2,4,i)
    cond=onOff>0&bothContr;    
    a=reshape(HC_pos(:,i,cond),700,sum(cond));
    [n,m]=max(a(200:400,:));
    a=reshape(HC_neg(:,i,cond),700,sum(cond));
    [n1,m1]=max(a(200:400,:));
    plot(m-m1,'.r')
    hold on
    
    cond=onOff<0&bothContr;    
    a=reshape(HC_pos(:,i,cond),700,sum(cond));
    [n,m]=min(a(200:400,:));
    a=reshape(HC_neg(:,i,cond),700,sum(cond));
    [n1,m1]=min(a(200:400,:));
    plot(m-m1,'.')
%     axis([1,200,-50,50])
end



figure
crit1=onOff>0;
crit2=ratio_area_HC;
crit3=zc_LC>50&zc_HC>50;
pl=wf;
baselined=1;
for i=1:8
    subplot(2,4,i)
    
    preC=crit1 & bothContr & crit3(i,:) & ~isnan(crit2(:,i))';
    
    cond=preC & (crit2(:,i)'<0.3);
    a=pl{i};
    if baselined
        a=a-repmat(mean(a(50:450,:)),4500,1);
    end
    plot(mean(a(:,cond),2),'r')  
    k=sum(cond);
    
    hold on
    
    cond=preC & (crit2(:,i)'>0.95);
    a=pl{i};
    if baselined
        a=a-repmat(mean(a(50:450,:)),4500,1);
    end
    plot(mean(a(:,cond),2),'b')    
    k1=sum(cond);
    
    cond=preC & (crit2(:,i)'<=0.95) & (crit2(:,i)'>=0.3);
    a=pl{i};
    if baselined
        a=a-repmat(mean(a(50:450,:)),4500,1);
    end
    plot(mean(a(:,cond),2),'g')
    k2=sum(cond);
    
    legend('Neg','Pos','Mixed')
    
    title([int2str(k),'  ',int2str(k1),'   ',int2str(k2),'  ', int2str(sum(preC))])
end





figure
for i=1:8
    subplot(2,4,i)
    cond=onOff>0&bothContr&zc_LC(i,:)>50&zc_HC(i,:)>50;
    plot(ratio_area_HC(cond,i),peak_pos_HC(cond,i)-peak_pos_LC(cond,i),'or')
    hold on
    cond=onOff<0&bothContr&zc_LC(i,:)>50&zc_HC(i,:)>50;
    plot(ratio_area_HC(cond,i),peak_pos_HC(cond,i)-peak_pos_LC(cond,i),'ob')
    axis([-1 1 -50 50])
end

figure
for i=1:8
    subplot(2,4,i)
    cond=onOff>0&bothContr&zc_LC(i,:)>50&zc_HC(i,:)>50;
    plot(ratio_area_HC(cond,i),peak_neg_HC(cond,i)-peak_neg_LC(cond,i),'or')
    hold on
    cond=onOff<0&bothContr&zc_LC(i,:)>50&zc_HC(i,:)>50;
    plot(ratio_area_HC(cond,i),peak_neg_HC(cond,i)-peak_neg_LC(cond,i),'ob')
    axis([-1 1 -50 50])
end

   

figure
for i=1:8
    subplot(2,4,i)
    cond=onOff>0&bothContr;
    plot(peak_neg_HC(cond,i),peak_neg_LC(cond,i),'or')
    hold on
    cond=onOff<0&bothContr;
    plot(peak_neg_HC(cond,i),peak_neg_LC(cond,i),'ob')
    
end



figure
for i=1:8
    subplot(2,4,i)
    cond=onOff>0&bothContr;
    plot(ratio_area_HC(cond,i),ratio_area_LC(cond,i),'or')
    hold on
    cond=onOff<0&bothContr;
    plot(ratio_area_HC(cond,i),ratio_area_LC(cond,i),'ob')
    axis([-1 1 -1 1])
end
