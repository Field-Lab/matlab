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


zc_HC_all=[];
zc_LC_all=[];
peak_HC_all=[];
peak_LC_all=[];

names_all=[];
onOff_all=[];
exp_codes=[];

bothContr=[];


for datesCNT=1:15
    date=dates{datesCNT}       
    
    path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];        
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_1',date])
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date],'names','onOff')
    names_all=[names_all; names'];
    onOff_all=[onOff_all onOff];
    exp_codes=[exp_codes ones(1,size(peak_HC,1))*datesCNT];
    
    if datesCNT<6
        bothContr=[bothContr zeros(1,size(peak_HC,1))];
    else
        bothContr=[bothContr ones(1,size(peak_HC,1))];
    end
    
    
    
    zc_HC_all=[zc_HC_all; zc_HC];
    zc_LC_all=[zc_LC_all; zc_LC];
    peak_HC_all=[peak_HC_all; peak_HC];
    peak_LC_all=[peak_LC_all; peak_LC];

end

zc_HC=zc_HC_all';
zc_LC=zc_LC_all';
peak_HC=peak_HC_all';
peak_LC=peak_LC_all';
names=names_all;
onOff=onOff_all;

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/summaryzc_late','peak_HC','peak_LC','bothContr','names','onOff','zc_HC','zc_LC','dates','exp_codes')




cd('/mnt/muench_data/user/alexandra/scripts')
clear
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summaryzc_late')

%% Figure M1

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
