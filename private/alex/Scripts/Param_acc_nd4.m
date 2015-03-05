clear

dates=cell(2,1);
dates{1}='20120902_1';
dates{2}='20120902_2';


zc_HC_all=[];
zc_LC_all=[];
peak_HC_all=[];
peak_LC_all=[];
wf=cell(3,1);
bf=cell(3,1);
frMean_HC_all=[];
frMean_LC_all=[];
frSTD_HC_all=[];
frSTD_LC_all=[];
frMean_spont_all=[];
frSTD_spont_all=[];
HC=[];
LC=[];

names_all=[];
onOff_all=[];
exp_codes=[];

bothContr=[];


for datesCNT=1:2
    date=dates{datesCNT}       
    
    path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];        
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])
    load([path2save,'quick'],'white_flash','black_flash','names')
    formula=formula(1:3);

    names_all=[names_all; names];
    
    exp_codes=[exp_codes ones(1,size(LinearFilter,4))*datesCNT];
    
    bothContr=[bothContr ones(1,size(LinearFilter,4))];
    
    zc_HC_all=[zc_HC_all; zc_HC(:,1:3)];
    zc_LC_all=[zc_LC_all; zc_LC(:,1:3)];
    peak_HC_all=[peak_HC_all; peak_HC(:,1:3)];
    peak_LC_all=[peak_LC_all; peak_LC(:,1:3)];

    frMean_HC_all=[frMean_HC_all frMean_HC(formula,:)];
    frMean_LC_all=[frMean_LC_all frMean_LC(formula,:)];
    
    frSTD_HC_all=[frSTD_HC_all frSTD_HC(formula,:)];
    frSTD_LC_all=[frSTD_LC_all frSTD_LC(formula,:)];

    frMean_spont_all=[frMean_spont_all frspont(formula,:)];
    frSTD_spont_all=[frSTD_spont_all frSTDspont(formula,:)];

    onOff_all=[onOff_all onOff];
    
    t=1;
    p=1;
    tt=size(HC,3)+1*~isempty(HC):size(LinearFilter,4)+size(HC,3)-1*isempty(HC);
    
    for i=4:5:15
        tmp=mean(white_flash(:,i:i+1,:),2);
        tmp=reshape(tmp,4500,length(onOff));
        wf{t}=[wf{t} tmp];
        
        tmp=mean(black_flash(:,i:i+1,:),2);
        tmp=reshape(tmp,4500,length(onOff));
        bf{t}=[bf{t} tmp];
        
        t=t+1;
    end
    for i=1:4:12
        HC(:,p,tt)=reshape(mean(mean(LinearFilter(:,1:2:end,i:i+3,:),2),3),500,size(LinearFilter,4));
        LC(:,p,tt)=reshape(mean(mean(LinearFilter(:,2:2:end,i:i+3,:),2),3),500,size(LinearFilter,4));
        p=p+1;
    end
    
end

frMean_HC=frMean_HC_all';
frMean_LC=frMean_LC_all';
frSTD_HC=frSTD_HC_all';
frSTD_LC=frSTD_LC_all';
frMean_spont=frMean_spont_all';
frSTD_spont=frSTD_spont_all';
zc_HC=zc_HC_all';
zc_LC=zc_LC_all';
peak_HC=peak_HC_all';
peak_LC=peak_LC_all';
names=names_all;
onOff=onOff_all;

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_nd4','peak_HC','peak_LC','HC','LC','bothContr','bf','wf','names','onOff','zc_HC','zc_LC','frMean_HC','frMean_LC','frSTD_HC','frSTD_LC','frMean_spont','frSTD_spont','dates','exp_codes')
clear all

%% Temporal development
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_nd4')
cd('/mnt/muench_data/user/alexandra/scripts')
clear

dates=cell(2,1);
dates{1}='20120902_1';
dates{2}='20120902_2';

HCon=[];
HCoff=[];
names_all=[];
onOff_all=[];
exp_codes=[];
for datesCNT=1:2
    date=dates{datesCNT}   
    
    
    path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
    load([path2save,'quick'],'white_flash','black_flash','names')
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])    
    
    names_all=[names_all names];    
    exp_codes=[exp_codes ones(1,size(black_flash,3))*datesCNT];    
    onOff_all=[onOff_all onOff];
    
    tt=size(HCon,3)+1*~isempty(HCon):sum(onOff>0)+size(HCon,3)-1*isempty(HCon);
    pp=size(HCoff,3)+1*~isempty(HCoff):sum(onOff<0)+size(HCoff,3)-1*isempty(HCoff);
    
    HCon(:,:,tt)=reshape(LinearFilter(:,1,1:12,onOff>0),500,12,sum(onOff>0));
    HCoff(:,:,pp)=reshape(LinearFilter(:,1,1:12,onOff<0),500,12,sum(onOff<0));
    
end
names=names_all;
onOff=onOff_all;
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/temporal_develop_LF_nd4.mat','HCon','HCoff','names','onOff','exp_codes')







%% Plot like M2
clear
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/temporal_develop_LF_nd4.mat','HCon','HCoff','names','onOff','exp_codes')

c=1;
for i=1:size(HCoff,3)
    HCoff(:,:,i)=HCoff(:,:,i)./repmat(sum(abs(HCoff(:,:,i))),500,1);
end
for i=1:3
    [~,peak_off(i,:)]=min(nanmean(HCoff(:,c:c+3,:),3));
    [~,peak_on(i,:)]=max(nanmean(HCon(:,c:c+3,:),3));
    for m=1:4
        zc_off(i,m)=find(nanmean(HCoff(mean(peak_off(i,:)):end,c+m-1,:),3)>0,1)+mean(peak_off(i,:))-1;
        zc_on(i,m)=find(nanmean(HCon(mean(peak_on(i,:)):end,c+m-1,:),3)<0,1)+mean(peak_on(i,:))-1;
    end
    c=c+4;
end

figure
c=5;
for i=1:3
    hold on
    plot(c:c+3,peak_on(i,:),'.r')
    plot(c:c+3,peak_off(i,:),'.b')
    plot(c:c+3,zc_on(i,:),'.m')
    plot(c:c+3,zc_off(i,:),'.c')
    c=c+7;        
end
legend({'peak ON','peak OFF','zero crossing ON','zero crossing OFF'})
axis([0 25 0 200])
set(gca,'xtick',6.5:7:25,'xticklabel',{'4','3','2'})
xlabel('ND')
ylabel('latency, ms')
title('Fig.M2. Latencies of peak and zero crossing for ON and OFF cells in temporal development. High contrast GWN')

