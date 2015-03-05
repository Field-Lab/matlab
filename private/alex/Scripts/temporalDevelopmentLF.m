cd('/mnt/muench_data/user/alexandra/scripts')
clear

dates=cell(6,1);
dates{1}='20130220';
dates{2}='20130220_1';
dates{3}='20130224';
dates{4}='20130225';
dates{5}='20130226';
dates{6}='20121023';
k=[];
for i=1:24:24*8
    k=[k i:2:i+17];
end
HCon=[];
HCoff=[];
names_all=[];
onOff_all=[];
exp_codes=[];
for datesCNT=1:6
    date=dates{datesCNT}    
    
    
    path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
    load([path2save,date,'_nonlinear'],'nonLinear')
    load([path2save,'quick'],'white_flash','black_flash','names')
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])    
    
    names_all=[names_all names];    
    exp_codes=[exp_codes ones(1,size(black_flash,3))*datesCNT];    
    onOff_all=[onOff_all onOff];
    
    tt=size(HCon,3)+1*~isempty(HCon):sum(onOff>0)+size(HCon,3)-1*isempty(HCon);
    pp=size(HCoff,3)+1*~isempty(HCoff):sum(onOff<0)+size(HCoff,3)-1*isempty(HCoff);

    switch size(LinearFilter,3)
        case 82           
            if strcmp(date,'20130220')
                HCon(:,:,tt)=reshape(LinearFilter(:,1,11:82,onOff>0),500,72,sum(onOff>0));
                HCoff(:,:,pp)=reshape(LinearFilter(:,1,11:82,onOff<0),500,72,sum(onOff<0));
            else
                HCon(:,:,tt)=reshape(LinearFilter(:,1,10:81,onOff>0),500,72,sum(onOff>0));
                HCoff(:,:,pp)=reshape(LinearFilter(:,1,10:81,onOff<0),500,72,sum(onOff<0));
            end

        case 342
            HCon(:,:,tt)=reshape(LinearFilter(:,1,k,onOff>0),500,72,sum(onOff>0));
            HCoff(:,:,pp)=reshape(LinearFilter(:,1,k,onOff<0),500,72,sum(onOff<0));
            
    end
    
end
names=names_all;
onOff=onOff_all;
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/temporal_develop_LF.mat','HCon','HCoff','names','onOff','exp_codes')


clear
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/temporal_develop_LF.mat','HCon','HCoff','names','onOff','exp_codes')


figure
c=1;
for i=1:size(HCoff,3)
    HCoff(:,:,i)=HCoff(:,:,i)./repmat(sum(abs(HCoff(:,:,i))),500,1);
end
for i=1:8
    subplot(2,4,i)
    hold on
    plot(mean(-HCoff(:,c:c+8,:),3),'b')
    [~,peak_off(i,:)]=min(nanmean(HCoff(:,c:c+8,:),3));
    for m=1:9
        zc_off(i,m)=find(nanmean(HCoff(mean(peak_off(i,:)):end,c+m-1,:),3)>0,1)+mean(peak_off(i,:))-1;
    end
    c=c+9;
end

c=1;
for i=1:size(HCon,3)
    HCon(:,:,i)=HCon(:,:,i)./repmat(sum(abs(HCon(:,:,i))),500,1);
end
for i=1:8
    subplot(2,4,i)
    plot(mean(HCon(:,c:c+8,:),3),'r')
    [~,peak_on(i,:)]=max(nanmean(HCon(:,c:c+8,:),3));
    for m=1:9
        zc_on(i,m)=find(nanmean(HCon(mean(peak_on(i,:)):end,c+m-1,:),3)<0,1)+mean(peak_on(i,:))-1;
    end
    c=c+9;
end



figure
c=1;
for i=1:8
    hold on
    plot(c:c+8,peak_on(i,:),'.r')
    plot(c:c+8,peak_off(i,:),'.b')
    plot(c:c+8,zc_on(i,:),'.m')
    plot(c:c+8,zc_off(i,:),'.c')
    c=c+12;
    axis([0 100 0 250])
end
