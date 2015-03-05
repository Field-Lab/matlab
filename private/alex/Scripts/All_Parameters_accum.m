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

gain_HC=[];
gain_LC=[];

zc_HC_all=[];
zc_LC_all=[];
peak_HC_all=[];
peak_LC_all=[];
wf=cell(8,1);
bf=cell(8,1);
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


for datesCNT=1:15
    date=dates{datesCNT}    
    
    
    path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];        
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])
    load([path2save,date,'_nonlinear'],'nonLinear')
    load([path2save,'quick'],'white_flash','black_flash','names')

    names_all=[names_all; names];
    
    exp_codes=[exp_codes ones(1,size(black_flash,3))*datesCNT];
    
    if size(black_flash,2)==72
        bothContr=[bothContr zeros(1,size(black_flash,3))*datesCNT];
    else
        bothContr=[bothContr ones(1,size(black_flash,3))];
    end
    
    zc_HC_all=[zc_HC_all; zc_HC];
    zc_LC_all=[zc_LC_all; zc_LC];
    peak_HC_all=[peak_HC_all; peak_HC];
    peak_LC_all=[peak_LC_all; peak_LC];

    frMean_HC_all=[frMean_HC_all frMean_HC(formula,:)];
    frMean_LC_all=[frMean_LC_all frMean_LC(formula,:)];
    
    frSTD_HC_all=[frSTD_HC_all frSTD_HC(formula,:)];
    frSTD_LC_all=[frSTD_LC_all frSTD_LC(formula,:)];

    frMean_spont_all=[frMean_spont_all frspont(formula,:)];
    frSTD_spont_all=[frSTD_spont_all frSTDspont(formula,:)];

    onOff_all=[onOff_all onOff];
    
    size(black_flash,2)
    t=1;
    
    p=1;
    tt=size(HC,3)+1*~isempty(HC):size(LinearFilter,4)+size(HC,3)-1*isempty(HC);
    
    switch size(black_flash,2)
        case 72
            for i=8:9:72
                tmp=mean(white_flash(:,i:i+1,:),2);
                tmp=reshape(tmp,4500,length(onOff));
                wf{t}=[wf{t} tmp];
                
                tmp=mean(black_flash(:,i:i+1,:),2);
                tmp=reshape(tmp,4500,length(onOff));
                bf{t}=[bf{t} tmp];
                
                t=t+1;
            end
            strt=10;
            if strcmp(date,'20130220')
                strt=11;
            end
            for i=strt:9:81
                HC(:,p,tt)=reshape(mean(LinearFilter(:,1,i:i+8,:),3),500,size(LinearFilter,4));
                LC(:,p,tt)=reshape(mean(LinearFilter(:,1,i:i+8,:),3),500,size(LinearFilter,4))*0;
                gain_HC(:,p,tt)=reshape(mean(nonLinear(:,1,i:i+8,:),3),6,size(LinearFilter,4));
                gain_LC(:,p,tt)=reshape(mean(nonLinear(:,1,i:i+8,:),3),6,size(LinearFilter,4))*0;
                p=p+1;
            end
        case 32
            for i=3:4:32
                tmp=mean(white_flash(:,i:i+1,:),2);
                tmp=reshape(tmp,4500,length(onOff));
                wf{t}=[wf{t} tmp];
                
                tmp=mean(black_flash(:,i:i+1,:),2);
                tmp=reshape(tmp,4500,length(onOff));
                bf{t}=[bf{t} tmp];
                
                t=t+1;
            end
            for i=5:4:36
                HC(:,p,tt)=reshape(mean(mean(LinearFilter(:,1:2:end,i:i+3,:),2),3),500,size(LinearFilter,4));
                LC(:,p,tt)=reshape(mean(mean(LinearFilter(:,2:2:end,i:i+3,:),2),3),500,size(LinearFilter,4));                
                gain_HC(:,p,tt)=reshape(mean(mean(nonLinear(:,1:2:end,i:i+3,:),2),3),6,size(LinearFilter,4));
                gain_LC(:,p,tt)=reshape(mean(mean(nonLinear(:,2:2:end,i:i+3,:),2),3),6,size(LinearFilter,4));
                p=p+1;
            end
        case 16
            for i=2:2:16
                tmp=white_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                wf{t}=[wf{t} tmp];
                
                tmp=black_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                bf{t}=[bf{t} tmp];
                
                t=t+1;
            end
            for i=1:2:16
                HC(:,p,tt)=reshape(mean(mean(LinearFilter(:,1:2:end,i:i+1,:),2),3),500,size(LinearFilter,4));
                LC(:,p,tt)=reshape(mean(mean(LinearFilter(:,2:2:end,i:i+1,:),2),3),500,size(LinearFilter,4));
                gain_HC(:,p,tt)=reshape(mean(mean(nonLinear(:,1:2:end,i:i+1,:),2),3),6,size(LinearFilter,4));
                gain_LC(:,p,tt)=reshape(mean(mean(nonLinear(:,2:2:end,i:i+1,:),2),3),6,size(LinearFilter,4));
                p=p+1;
            end
        case 17
            for i=2:2:16
                tmp=white_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                wf{t}=[wf{t} tmp];
                
                tmp=black_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                bf{t}=[bf{t} tmp];

                t=t+1;
            end
            for i=1:2:16
                HC(:,p,tt)=reshape(mean(mean(LinearFilter(:,1:2:end,i:i+1,:),2),3),500,size(LinearFilter,4));
                LC(:,p,tt)=reshape(mean(mean(LinearFilter(:,2:2:end,i:i+1,:),2),3),500,size(LinearFilter,4));
                gain_HC(:,p,tt)=reshape(mean(mean(nonLinear(:,1:2:end,i:i+1,:),2),3),6,size(LinearFilter,4));
                gain_LC(:,p,tt)=reshape(mean(mean(nonLinear(:,2:2:end,i:i+1,:),2),3),6,size(LinearFilter,4));
                p=p+1;
            end
        case 24
            for i=3:3:24
                tmp=white_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                wf{t}=[wf{t} tmp];
                
                tmp=black_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                bf{t}=[bf{t} tmp];

                t=t+1;
            end
            for i=1:3:24
                HC(:,p,tt)=reshape(mean(mean(LinearFilter(:,1:2:end,i:i+2,:),2),3),500,size(LinearFilter,4));
                LC(:,p,tt)=reshape(mean(mean(LinearFilter(:,2:2:end,i:i+2,:),2),3),500,size(LinearFilter,4));
                gain_HC(:,p,tt)=reshape(mean(mean(nonLinear(:,1:2:end,i:i+2,:),2),3),6,size(LinearFilter,4));
                gain_LC(:,p,tt)=reshape(mean(mean(nonLinear(:,2:2:end,i:i+2,:),2),3),6,size(LinearFilter,4));
                p=p+1;
            end
        case 56
            for i=4:5:40
                tmp=mean(white_flash(:,i:i+1,:),2);
                tmp=reshape(tmp,4500,length(onOff));
                wf{t}=[wf{t} tmp];
                
                tmp=mean(black_flash(:,i:i+1,:),2);
                tmp=reshape(tmp,4500,length(onOff));
                bf{t}=[bf{t} tmp];
                
                t=t+1;
            end
            for i=1:4:32
                HC(:,p,tt)=reshape(mean(mean(LinearFilter(:,1:2:end,i:i+3,:),2),3),500,size(LinearFilter,4));
                LC(:,p,tt)=reshape(mean(mean(LinearFilter(:,2:2:end,i:i+3,:),2),3),500,size(LinearFilter,4));
                gain_HC(:,p,tt)=reshape(mean(mean(nonLinear(:,1:2:end,i:i+3,:),2),3),6,size(LinearFilter,4));
                gain_LC(:,p,tt)=reshape(mean(mean(nonLinear(:,2:2:end,i:i+3,:),2),3),6,size(LinearFilter,4));
                p=p+1;
            end
        case 8
            for i=1:8
                tmp=white_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                wf{t}=[wf{t} tmp];
                
                tmp=black_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                bf{t}=[bf{t} tmp];
                
                t=t+1;
            end
            for i=1:6:6*8
                HC(:,p,tt)=reshape(mean(mean(LinearFilter(:,1:2:end,i:i+5,:),2),3),500,size(LinearFilter,4));
                LC(:,p,tt)=reshape(mean(mean(LinearFilter(:,2:2:end,i:i+5,:),2),3),500,size(LinearFilter,4));
                gain_HC(:,p,tt)=reshape(mean(mean(nonLinear(:,1:2:end,i:i+5,:),2),3),6,size(LinearFilter,4));
                gain_LC(:,p,tt)=reshape(mean(mean(nonLinear(:,2:2:end,i:i+5,:),2),3),6,size(LinearFilter,4));
                p=p+1;
            end
        case 7
            t=2;
            for i=1:7
                tmp=white_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                wf{t}=[wf{t} tmp];
                
                tmp=black_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                bf{t}=[bf{t} tmp];
                
                t=t+1;
            end
            wf{1}=[wf{1} zeros(4500,length(onOff))];
            bf{1}=[bf{1} zeros(4500,length(onOff))];
            p=2;
            for i=1:6:6*7
                HC(:,p,tt)=reshape(mean(mean(LinearFilter(:,1:2:end,i:i+5,:),2),3),500,size(LinearFilter,4));
                LC(:,p,tt)=reshape(mean(mean(LinearFilter(:,2:2:end,i:i+5,:),2),3),500,size(LinearFilter,4));
                gain_HC(:,p,tt)=reshape(mean(mean(nonLinear(:,1:2:end,i:i+5,:),2),3),6,size(LinearFilter,4));
                gain_LC(:,p,tt)=reshape(mean(mean(nonLinear(:,2:2:end,i:i+5,:),2),3),6,size(LinearFilter,4));
                p=p+1;
            end
        case 14
            for i=1:8
                tmp=white_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                wf{t}=[wf{t} tmp];
                
                tmp=black_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                bf{t}=[bf{t} tmp];

                t=t+1;
            end
            for i=1:24:24*8
                HC(:,p,tt)=reshape(mean(LinearFilter(:,1,i:2:i+23,:),3),500,size(LinearFilter,4));
                LC(:,p,tt)=reshape(mean(LinearFilter(:,1,(i+1):2:i+23,:),3),500,size(LinearFilter,4));
                gain_HC(:,p,tt)=reshape(mean(nonLinear(:,1,i:2:i+23,:),3),6,size(LinearFilter,4));
                gain_LC(:,p,tt)=reshape(mean(nonLinear(:,1,(i+1):2:i+23,:),3),6,size(LinearFilter,4));
                p=p+1;
            end
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

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all','gain_HC','gain_LC','peak_HC','peak_LC','HC','LC','bothContr','bf','wf','names','onOff','zc_HC','zc_LC','frMean_HC','frMean_LC','frSTD_HC','frSTD_LC','frMean_spont','frSTD_spont','dates','exp_codes')
clear all


%% Plots
clear
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_all')


% Plot every cell separately
path2save='/mnt/muench_data/data/alexandra/MEA_data/analysis/Quick_0_cells/';
if ~exist(path2save,'dir')
    mkdir(path2save)
end
figure
set(gcf,'position',[20 50        1639         900])
nds='87654321'
for cnt=1:size(bf{1},2)
    if onOff(cnt)==0
        maxFR=[];
        
        for i=1:8
            subplot(4,2,i)
            hold off
            plot(bf{i}(:,cnt),'r','LineWidth',2)
            hold on
            maxFR=[maxFR max(bf{i}(:,cnt))];
            plot(4701:9200,wf{i}(:,cnt),'LineWidth',2)
            maxFR=[maxFR max(wf{i}(:,cnt))];
        end
        k=max(maxFR)*1.1;
        for i=1:8
            subplot(4,2,i)
            axis([0 9200 0 k]);
            line([500,500],[0 k],'color','k')
            line([2500,2500],[0 k],'color','k')
            line([5200,5200],[0 k],'color','k')
            line([7200,7200],[0 k],'color','k')
            title(['ND',nds(i)])
        end
        subplot('Position',[0.5 0.96 0.00001 0.00001])
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        title([int2str(cnt),'   ',names{cnt},'   RED - BLACK FLASH, BLUE - WHITE FLASH'],'FontSize',12,'FontWeight','bold','Interpreter','None')
        saveas(gcf,[path2save,names{cnt},'.png'])
    end
end







figure
d=find(onOff<0);
% m=1:200;
for i=1:8
    subplot(4,2,i)
    exp_codes(d(m));
    imagesc([bf{i}(:,d(m))' wf{i}(:,d(m))' zeros(200,1) zeros(200,1)+200 repmat(exp_codes(d(m))',1,500)*13])
    max(max([bf{i}(:,d(m))' wf{i}(:,d(m))']));
    title(['ND',nds(i)])
end    


k=wf{4}(:,onOff<0)';
k=k-repmat(mean(k(:,50:450),2),1,4500);
figure
plot(k')
[~, m]=sort(mean(k(:,1100:1600)'));

i=3
k=bf{i}(:,onOff<0)';
k=k-repmat(mean(k(:,50:450),2),1,4500);
[~, m]=sort(mean(k(:,3100:3500)'));
kwf=wf{i}(:,onOff<0)';
kwf=kwf-repmat(mean(kwf(:,50:450),2),1,4500);
figure
imagesc([k(m,:) kwf(m,:) zeros(200,1) zeros(200,1)+200 repmat(exp_codes(d(m))',1,500)*13])
figure
for i=1:8
    subplot(4,2,i)
    k=bf{i}(:,onOff<0)';
    k=k-repmat(mean(k(:,50:450),2),1,4500);
    kwf=wf{i}(:,onOff<0)';
    kwf=kwf-repmat(mean(kwf(:,50:450),2),1,4500);
    plot([k(m(end-30:end),:) kwf(m(end-30:end),:)]')
end

i=4
k=bf{i}(:,onOff<0)';
k=k-repmat(mean(k(:,50:450),2),1,4500);
kwf=wf{i}(:,onOff<0)';
kwf=kwf-repmat(mean(kwf(:,50:450),2),1,4500);
k=[k kwf];
a=corr(k');
a(a==1)=-1;
d=max(a(:));
[r,c]=find(a==d,1);

sorted(1)=r;
sorted(2)=c;
p(1)=d;
for j=2:199
    a=corr(k');
    a(a==1)=-1;
    a(:,sorted(1:j-1))=-1;
    a(sorted(1:j-1),:)=-1;
    [d,r]=max(a(:,sorted(j)));
    p(j)=d;
    sorted(j+1)=r;
end
m=sorted;
i=4
k=bf{i}(:,onOff<0)';
k=k-repmat(mean(k(:,50:450),2),1,4500);
kwf=wf{i}(:,onOff<0)';
kwf=kwf-repmat(mean(kwf(:,50:450),2),1,4500);
figure
d=find(onOff<0);
imagesc([k(m,:) kwf(m,:) zeros(200,1) zeros(200,1)+200 repmat(exp_codes(d(m))',1,500)*13])



for i=1:8
    k=bf{i}(:,onOff<0)';
    k=k-repmat(mean(k(:,50:450),2),1,4500);
    kwf=wf{i}(:,onOff<0)';
    kwf=kwf-repmat(mean(kwf(:,50:450),2),1,4500);
    
    offBF(i,:)=max(k(:,500:750)')-max(k(:,50:350)');
    offWF(i,:)=max(kwf(:,2500:2750)')-max(kwf(:,50:350)');
    rebBF(i,:)=max(k(:,3050:3500)')-max(k(:,50:350)');
    rebWF(i,:)=max(kwf(:,1050:1500)')-max(kwf(:,50:350)');
    onBF(i,:)=max(k(:,2600:2750)')-max(k(:,50:350)');
    onWF(i,:)=max(kwf(:,600:750)')-max(kwf(:,50:350)');
end
    
figure
[a,b]=sort(rebBF(4,:));
imagesc([rebBF(2:7,b)' zeros(200,1)-40 rebWF(2:7,b)' zeros(200,1)-40 onBF(2:7,b)' zeros(200,1)-40 onWF(2:7,b)'...
    zeros(200,1)-40 offBF(2:7,b)' zeros(200,1)-40 offWF(2:7,b)'])

   imagesc([rebBF(2:7,b)'>2 zeros(200,1)-0.5 rebWF(2:7,b)'>2 zeros(200,1)-0.5 onBF(2:7,b)'>2 zeros(200,1)-0.5 ...
       onWF(2:7,b)'>2 zeros(200,1)-0.5 offBF(2:7,b)'>2 zeros(200,1)-0.5 offWF(2:7,b)'>2])

   a=[rebBF(3:6,b)'>2 rebWF(3:6,b)'>2 onBF(3:6,b)'>2 onWF(3:6,b)'>2];

   
%% OFF cells


%% Figs Z

% ON and OFF population mean responses for fig Z and R
figure
nds='87654321'
ind_start1= 700;%wf, OFF cells: 700 for ON-OFF, 1200 for rebound, 2500 for OFF response
ind_start2= 2700; % bf, OFF cells: 700 for OFF response, 2700 for ON-OFF, 3200 for rebound
for i=1:8
    subplot(4,2,i)
    tmp=wf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    hold off
    plot(mean(tmp,2),'color','k','linewidth',2)
    %indices
%     [a,b]=max(tmp(ind_start1-100:ind_start1+100,:));
%     ampl1(i,1)=mean(a);
%     ampl1(i,2)=std(a)/sqrt(length(a));
% 
%     lat1(i,1)=mean(b);
%     lat1(i,2)=std(b)/sqrt(length(b));
%     
%     
    hold on
    tmp=bf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    plot(4701:9200,mean(tmp,2),'color','k','linewidth',2)

    %indices
%     [a,b]=max(tmp(ind_start2-100:ind_start2+100,:));
%     ampl2(i,1)=mean(a);
%     ampl2(i,2)=std(a)/sqrt(length(a));
% 
%     lat2(i,1)=mean(b);
%     lat2(i,2)=std(b)/sqrt(length(b));    
%     
    
    
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


% 
%     axis([1 9200 -20 50])
%     line([1 4500],[0,0],'color','k')
%     line([4700 9200],[0,0],'color','k')
%     
%     line([1 500],[40,40],'color','k')
%     line([500 2500],[45,45],'color','k')
%     line([2500 4500],[40,40],'color','k')
%     line([500 500],[40,45],'color','k')
%     line([2500 2500],[40,45],'color','k')
%     
%     line([1 500]+4700,[40,40],'color','k')
%     line([500 2500]+4700,[35,35],'color','k')
%     line([2500 4500]+4700,[40,40],'color','k')
%     line([500 500]+4700,[35,40],'color','k')
%     line([2500 2500]+4700,[35,40],'color','k')
    title(['ND',nds(i)])
end




load('/mnt/muench_data/data/alexandra/MEA_data/analysis/OFFcellsClust.txt')
manual=zeros(520,1)-1;
rebSt=0;onSt=0;
on=[0 0 0];reb=[0 0 0];
rebStIs=0;onStIs=0;
for i=1:200
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
    end
    if a(2)=='1'
        on(2)=on(2)+1;
    end
    if a(3)=='1'
        on(3)=on(3)+1;
    end
    if a(4)=='1'
        reb(1)=reb(1)+1;
    end
    if a(5)=='1'
        reb(2)=reb(2)+1;
    end
    if a(6)=='1'
        reb(3)=reb(3)+1;
    end
    manual(OFFcellsClust(i,1))=bin2dec(a);    
end

figure 
hist(manual,64)
  
figure
clear means
a=manual;
cc=1;tt=1;
for code=0:63
    sum(a==code);
    if sum(a==code)>12
        for i=2:7
            if i==2
                code
                exx=exp_codes(a==code)
            end
            subplot(6,9,cc)
            k=bf{i}';
            k=k-repmat(mean(k(:,50:450),2),1,4500);
            kwf=wf{i}';
            kwf=kwf-repmat(mean(kwf(:,50:450),2),1,4500);
            plot([k(a==code,:) kwf(a==code,:)]')
            axis tight
            hold on
            means(1+4500*(i-2):4500+4500*(i-2),tt)=mean(k(a==code,:));
            plot(mean([k(a==code,:) kwf(a==code,:)]),'color','k','linewidth',2)
            title([dec2bin(code,6),'   ',int2str(sum(a==code)),' ',int2str(length(unique(exx)))])
            
            
            cc=cc+1;
        end
        tt=tt+1;
    end
end



a=manual;
a_list=[];nums=zeros(64,1);
for code=0:63
    nums(code+1)=sum(a==code);
    a_list=[a_list; dec2bin(code,6)];
end

%% ON cells

% R1 Mean and std of the response to white flash, ON cells
figure
for i=3:5
    subplot(1,3,i-2)
    tmp=wf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    b=mean(tmp');
    a=std(tmp');
    plot(b,'k','Linewidth',2)
    hold on
    plot(b-a,'k','Linewidth',1)
    plot(b+a,'k','Linewidth',1)
    title(['ND', nds(i)])
    axis([0 4500 -30 85])
    line([1 4500],[0,0],'color','k')
   
    line([1 500],[70,70],'color','k')
    line([500 2500],[75,75],'color','k')
    line([2500 4500],[70,70],'color','k')
    line([500 500],[75,70],'color','k')
    line([2500 2500],[75,70],'color','k')
    
end



figure
for i=3:5
    subplot(3,1,i-2)    
    k=wf{i}(:,onOff>0)';
    p=mean(k(:,50:450),2);
    p1=mean(kbf(:,600:1000),2);
    k=k-repmat(mean(k(:,50:450),2),1,4500);    
    t=mean(k(:,2590:2690),2)<-3;
    kbf=bf{i}(:,onOff>0)';
    kbf=kbf-repmat(mean(kbf(:,50:450),2),1,4500);
    t1=mean(kbf(:,600:1000),2)<-3;
    sum(t1&t)
    plot([k(t&t1,:) kbf(t&t1,:)]')    
    title([int2str(sum(t&t1)),' cells, sFR mean ', num2str(mean(p(t&t1))),' min ', num2str(min(p(t&t1)))])


%     axis([4750 7000 -50 100])
%     axis([2250 4500 0 150])
end


clear offBF offWF rebBF rebWF onBF onWF endofPrefBF
for i=1:8
    k=wf{i}(:,onOff>0)';
    k=k-repmat(mean(k(:,50:450),2),1,4500);
    kbf=bf{i}(:,onOff>0)';
    spont(i,:)=mean(kbf(:,50:450)');
    kbf=kbf-repmat(mean(kbf(:,50:450),2),1,4500);
    
    inhWF(i,:)=mean(k(:,2550:3000)');
    rebWF(i,:)=mean(k(:,3000:3500)');
    lateWF(i,:)=mean(k(:,3500:4500)');
    inhBF(i,:)=mean(kbf(:,500:1000)');
    rebBF(i,:)=mean(kbf(:,1000:1500)');
    lateBF(i,:)=mean(kbf(:,1500:2450)');
end
figure
for i=3:5
    subplot(1,3,i-2)
    plot(spont(i,:),lateBF(i,:),'o')
    hold on
    plot(spont(i,:),lateWF(i,:),'or')
end

figure
for i=1:8
    subplot(2,4,i)    
    k=wf{i}(:,onOff>0)';
%     k=k-repmat(mean(k(:,50:450),2),1,4500);
    kbf=bf{i}(:,onOff>0)';
%     kbf=kbf-repmat(mean(kbf(:,50:450),2),1,4500);
    plot([k kbf]')
    hold on 
    plot(mean([k kbf]),'k','linewidth',2)
    axis([4750 7000 0 100])
%     axis([2250 4500 0 150])
end
figure
col='rgb'
for i=3:5
    k=wf{i}(:,onOff>0)';
    plot(mean(k(:,50:450),2),'.','color',col(i-2))
    hold on
end

figure
for i=3:5
    subplot(3,3,i-2)    
    k=wf{i}(:,onOff>0)';
    p=mean(k(:,50:450),2);
    k=k-repmat(mean(k(:,50:450),2),1,4500);    
    t=mean(k(:,2590:2690),2)<-10;
    kbf=bf{i}(:,onOff>0)';
    kbf=kbf-repmat(mean(kbf(:,50:450),2),1,4500);
    t1=mean(kbf(:,600:1000),2)<-5;
    sum(t1&t)
    plot([k(t1&t,:) kbf(t1&t,:)]')
    hold on 
    plot(mean([k(t1&t,:) kbf(t1&t,:)]),'k','linewidth',2)
    
    line([500 500],[-50,100],'color','k')
    line([2500 2500],[-50,100],'color','k')
    line([5000 5000],[-50,100],'color','k')
    line([7000 7000],[-50,100],'color','k')
    line([0 9000],[0,0],'color','k')
    title([num2str(sum(t&t1&(p>5))/sum(p>5)*100),'% of ',int2str(sum(p>5)),' cells with spont FR>5 or ',num2str(sum(t&t1)/249*100),'% of 249 cells'])

%     axis([4750 7000 -50 100])
%     axis([2250 4500 0 150])
end
for i=3:5
    subplot(3,3,i+1)    
    k=wf{i}(:,onOff>0)';
    p=mean(k(:,50:450),2);
    k=k-repmat(mean(k(:,50:450),2),1,4500);    
    t=mean(k(:,2590:2690),2)>=-10&mean(k(:,2590:2690),2)<2;
    kbf=bf{i}(:,onOff>0)';
    kbf=kbf-repmat(mean(kbf(:,50:450),2),1,4500);
    t1=mean(kbf(:,600:1000),2)>=-5&mean(kbf(:,600:1000),2)<2;
    sum(t1&t)
    plot([k(t1&t,:) kbf(t1&t,:)]')
    hold on 
    plot(mean([k(t1&t,:) kbf(t1&t,:)]),'k','linewidth',2)
    
    line([500 500],[-50,100],'color','k')
    line([2500 2500],[-50,100],'color','k')
    line([5000 5000],[-50,100],'color','k')
    line([7000 7000],[-50,100],'color','k')
    line([0 9000],[0,0],'color','k')
    title([num2str(sum(t&t1&(p>5))/sum(p>5)*100),'% of ',int2str(sum(p>5)),' cells with spont FR>5 or ',num2str(sum(t&t1)/249*100),'% of 249 cells'])


%     axis([4750 7000 -50 100])
%     axis([2250 4500 0 150])
end
for i=3:5
    subplot(3,3,i+4)    
    k=wf{i}(:,onOff>0)';
    p=mean(k(:,50:450),2);
    k=k-repmat(mean(k(:,50:450),2),1,4500);    
    t=mean(k(:,2590:2690),2)>=2;
    kbf=bf{i}(:,onOff>0)';
    kbf=kbf-repmat(mean(kbf(:,50:450),2),1,4500);
    t1=mean(kbf(:,600:1000),2)>=2;
    sum(t1&t)
    plot([k(t1&t,:) kbf(t1&t,:)]')
    hold on 
    plot(mean([k(t1&t,:) kbf(t1&t,:)]),'k','linewidth',2)
    
    line([500 500],[-50,100],'color','k')
    line([2500 2500],[-50,100],'color','k')
    line([5000 5000],[-50,100],'color','k')
    line([7000 7000],[-50,100],'color','k')
    line([0 9000],[0,0],'color','k')
    title([num2str(sum(t&t1&(p>5))/sum(p<2)*100),'% of ',int2str(sum(p<2)),' cells with spont FR<2 or ',num2str(sum(t&t1)/249*100),'% of 249 cells'])


%     axis([4750 7000 -50 100])
%     axis([2250 4500 0 150])
end


figure
for i=3:5
    subplot(3,1,i-2)    
    k=wf{i}(:,onOff>0)';
    p=mean(k(:,50:450),2);
    p1=mean(kbf(:,600:1000),2);
    k=k-repmat(mean(k(:,50:450),2),1,4500);    
    t=mean(k(:,2590:2690),2)<-3;
    kbf=bf{i}(:,onOff>0)';
    kbf=kbf-repmat(mean(kbf(:,50:450),2),1,4500);
    t1=mean(kbf(:,600:1000),2)<-3;
    sum(t1&t)
    imagesc([k(t&t1,:) kbf(t&t1,:)])    
    title([int2str(sum(t&t1)),' cells, sFR mean ', num2str(mean(p(t&t1))),' min ', num2str(min(p(t&t1)))])


%     axis([4750 7000 -50 100])
%     axis([2250 4500 0 150])
end
figure
for i=1:7
    subplot(2,7,i)    
    k=wf{i}(:,onOff>0)';
    p=mean(k(:,50:450),2);
    k=k-repmat(mean(k(:,50:450),2),1,4500);
    p1=mean(k(:,2590:2800),2);
    p1late=mean(k(:,3150:3400),2);
    t=mean(k(:,2590:2800),2)<-3;
    
    
    kbf=bf{i}(:,onOff>0)';
    pw=mean(kbf(:,50:450),2);
    kbf=kbf-repmat(mean(kbf(:,50:450),2),1,4500);
    p2=mean(kbf(:,590:800),2);
    p2late=mean(kbf(:,1150:1400),2);
    t1=mean(kbf(:,590:800),2)<-3;
    sum(t1&t)
    plot(p,p1,'o') 
    hold on
    plot(p,p1late,'or')
    line([0,50],[0,-50],'color','k')
%     hold on
%     plot(p(t&t1),p1(t&t1),'.r')  
    subplot(2,7,i+7)
    plot(pw,p2,'o')
    hold on
    plot(pw,p2late,'or')
    line([0,50],[0,-50],'color','k')
%     hold on
%     plot(p(t&t1),p2(t&t1),'.r')      
    title([int2str(sum(t&t1)),' cells, sFR mean ', num2str(mean(p(t&t1))),' min ', num2str(min(p(t&t1)))])


%     axis([4750 7000 -50 100])
%     axis([2250 4500 0 150])
end



a=[rebBF(3:5,:)'>bord onBF(3:5,:)'>bord];
imagesc(a)
a=int2str(a);
a=a(:,1:3:end);
a=bin2dec(a);
%    hist(a,max(a))
%    plot(sort(a),'.')
[~,ind]=sort(a);
a=[rebBF(3:5,ind)'>bord(ind,:) onBF(3:5,ind)'>bord(ind,:)];
figure
imagesc(a)



















clear k
for i=2:7
    k(i-1,1)=corr(rebBF(i,:)', rebWF(i,:)');
    k(i-1,2)=corr(onBF(i,:)', onWF(i,:)');
    k(i-1,3)=corr(offBF(i,:)', offWF(i,:)');
    k(i-1,4)=corr(rebBF(i,:)', onBF(i,:)');
    k(i-1,5)=corr(rebWF(i,:)', onWF(i,:)');
    k(i-1,6)=corr(rebBF(i,:)', rebBF(i+1,:)');
end
k=[(7:-1:2)' k];

imagesc((offBF(2:7,:)-onBF(2:7,:))'./(offBF(2:7,:)+onBF(2:7,:))')
figure
for i=1:8
    subplot(3,3,i)
    plot((offBF(i,:)-onBF(i,:))'./(offBF(i,:)+onBF(i,:))')
end




[a,b]=sort(rebWF(4,:));
imagesc(rebWF(:,b)')

[a,b]=sort(onBF(4,:));
imagesc(onBF(:,b)')

[a,b]=sort(onWF(4,:));
imagesc(onWF(:,b)')
    
    
    

k=onOff<0;

for i=1:8
    sp(i,:)=mean(bf{i}(50:450,k));
    spwf(i,:)=mean(wf{i}(50:450,k));
    beg(i,:)=max(bf{i}(500:1000,k));
    fin(i,:)=max(bf{i}(2000:2500,k));
    postbeg(i,:)=max(bf{i}(2500:3000,k));
    postfin(i,:)=max(bf{i}(4000:4500,k));
end
[~,m]=sort(beg(3,:));
plot(sp(2:6,m)'-beg(2:6,m)')

figure
nds='87654321'
nd1=[2 3 3 4];
nd2=[3 4 5 5];
range=2501:3000;
for l=1:4
    for j=1:520
        if onOff(j)<0
            p=mean(reshape(bf{nd1(l)}(range,j),25,length(range)/25));
            p1=mean(reshape(bf{nd2(l)}(range,j),25,length(range)/25));
            a(j)=corr(p',p1');
        elseif onOff(j)>0
            p=mean(reshape(wf{nd1(l)}(range,j),25,length(range)/25));
            p1=mean(reshape(wf{nd2(l)}(range,j),25,length(range)/25));
            a(j)=corr(p',p1');
        end
    end
    subplot(2,2,l)
    plot(sort(a(onOff>0)),'r')
    hold on
    plot(sort(a(onOff<0)))
    line([0 250],[0.5, 0.5],'color','k')
    title(['ND',nds(nd1(l)),' and ND',nds(nd2(l))])
end