cd('/mnt/muench_data/user/alexandra/scripts')

clear

dates=cell(12,1);
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


wf=cell(8,1);
bf=cell(8,1);

names_all=[];
onOff_all=[];
exp_codes=[];

for datesCNT=1:12
    date=dates{datesCNT}   

    path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];        
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])
    load([path2save,'quick'],'white_flash','black_flash','names')

    names_all=[names_all; names];
    
    exp_codes=[exp_codes ones(1,size(black_flash,3))*datesCNT];
    
    
    onOff_all=[onOff_all onOff];
    
    size(black_flash,2)
    t=1;
    
    p=1;
    
    switch size(black_flash,2)
        case 72
            for i=1:9:72
                tmp=white_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                wf{t}=[wf{t} tmp];
                
                tmp=black_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                bf{t}=[bf{t} tmp];
                
                t=t+1;
            end
        case 32
            for i=1:4:32
                tmp=mean(white_flash(:,i,:),2);
                tmp=reshape(tmp,4500,length(onOff));
                wf{t}=[wf{t} tmp];
                
                tmp=mean(black_flash(:,i,:),2);
                tmp=reshape(tmp,4500,length(onOff));
                bf{t}=[bf{t} tmp];
                
                t=t+1;
            end

        case 16
            for i=1:2:16
                tmp=white_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                wf{t}=[wf{t} tmp];
                
                tmp=black_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                bf{t}=[bf{t} tmp];
                
                t=t+1;
            end
        case 17
            for i=1:2:16
                tmp=white_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                wf{t}=[wf{t} tmp];
                
                tmp=black_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                bf{t}=[bf{t} tmp];

                t=t+1;
            end
        case 24
            for i=1:3:24
                tmp=white_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                wf{t}=[wf{t} tmp];
                
                tmp=black_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                bf{t}=[bf{t} tmp];

                t=t+1;
            end
        case 56
            for i=1:5:40
                tmp=mean(white_flash(:,i,:),2);
                tmp=reshape(tmp,4500,length(onOff));
                wf{t}=[wf{t} tmp];
                
                tmp=mean(black_flash(:,i,:),2);
                tmp=reshape(tmp,4500,length(onOff));
                bf{t}=[bf{t} tmp];                
                t=t+1;
            end
    end
    
end

names=names_all;
onOff=onOff_all;

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_quick_firsts','bf','wf','names','onOff','dates','exp_codes')

clear all

dates=cell(12,1);
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


wf=cell(8,1);
bf=cell(8,1);

names_all=[];
onOff_all=[];
exp_codes=[];

for datesCNT=1:12
    date=dates{datesCNT}   

    path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];        
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])
    load([path2save,'quick'],'white_flash','black_flash','names')

    names_all=[names_all; names];
    
    exp_codes=[exp_codes ones(1,size(black_flash,3))*datesCNT];
    
    
    onOff_all=[onOff_all onOff];
    
    size(black_flash,2)
    t=1;
    
    p=1;
    
    switch size(black_flash,2)
        case 72
            for i=9:9:72
                tmp=white_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                wf{t}=[wf{t} tmp];
                
                tmp=black_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                bf{t}=[bf{t} tmp];
                
                t=t+1;
            end
        case 32
            for i=4:4:32
                tmp=white_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                wf{t}=[wf{t} tmp];
                
                tmp=black_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                bf{t}=[bf{t} tmp];
                
                t=t+1;
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
        case 56
            for i=5:5:40
                tmp=white_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                wf{t}=[wf{t} tmp];
                
                tmp=black_flash(:,i,:);
                tmp=reshape(tmp,4500,length(onOff));
                bf{t}=[bf{t} tmp];                
                t=t+1;
            end
    end
    
end

names=names_all;
onOff=onOff_all;

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_quick_lasts','bf','wf','names','onOff','dates','exp_codes')
clear all

%% Plot
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_quick_firsts')
bf1=bf;
wf1=wf;
load('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_quick_lasts')

dataOfffirst=zeros(4500,2,8);
dataOfflast=zeros(4500,2,8);
for i=1:8
    subplot(4,2,i)
    hold on
    
    tmp=wf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    plot(mean(tmp,2),'color','b','linewidth',2)
    dataOfflast(:,1,i)=mean(tmp,2);
    k=tmp;
        
    tmp=bf{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    plot(4701:9200,mean(tmp,2),'color','b','linewidth',2)
    dataOfflast(:,2,i)=mean(tmp,2);

        
        
    tmp=wf1{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    plot(mean(tmp,2),'color','r','linewidth',2)
    dataOfffirst(:,1,i)=mean(tmp,2);
    k=tmp;
    
    tmp=bf1{i}(:,onOff<0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    plot(4701:9200,mean(tmp,2),'color','r','linewidth',2)
    dataOfffirst(:,2,i)=mean(tmp,2);
    
    a(1:183,i)=max(k(900:1500,:));
    b(1:183,i)=max(k1(900:1500,:));
    for j=1:183
        pp(j,i)=corr(k(900:1500,j),k1(900:1500,j));
    end

end

dataOnfirst=zeros(4500,2,8);
dataOnlast=zeros(4500,2,8);
figure
for i=1:8
    subplot(4,2,i)
    hold on
    
    tmp=wf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    plot(mean(tmp,2),'color','b','linewidth',2)
    dataOnlast(:,1,i)=mean(tmp,2);
    tmp=bf{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    plot(4701:9200,mean(tmp,2),'color','b','linewidth',2)
    dataOnlast(:,2,i)=mean(tmp,2);

    
    tmp=wf1{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    plot(mean(tmp,2),'color','r','linewidth',2)
    dataOnfirst(:,1,i)=mean(tmp,2);
    tmp=bf1{i}(:,onOff>0);
    tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
    plot(4701:9200,mean(tmp,2),'color','r','linewidth',2)
    dataOnfirst(:,2,i)=mean(tmp,2);
    
    
end
    
save('/mnt/muench_data/data/alexandra/MEA_data/analysis/data_for_Q','dataOfffirst','dataOfflast','dataOnfirst','dataOnlast')


nanmean(pp)
nanstd(pp)

nanmean(a)
nanmean(b)
for i=1:8
    c(i)=ranksum(a(:,i),b(:,i));
end






clear all

dates=cell(5,1);
dates{1}='20130220';
dates{2}='20130220_1';
dates{3}='20130224';
dates{4}='20130225';
dates{5}='20130226';


wf=[];
bf=[];

names_all=[];
onOff_all=[];
l=1;

for datesCNT=1:5
    date=dates{datesCNT}   

    path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];        
    load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date],'onOff')
    load([path2save,'quick'],'white_flash','black_flash','names')

    names_all=[names_all; names];
    onOff_all=[onOff_all onOff];
    for i=1:72
        wf(:,l:l+size(onOff,2)-1,i)=reshape(white_flash(:,i,:),4500,size(onOff,2));
        bf(:,l:l+size(onOff,2)-1,i)=reshape(black_flash(:,i,:),4500,size(onOff,2));
    end
    l=l+size(onOff,2);
    
end

names=names_all;
onOff=onOff_all;

save('/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_quick_temporal','bf','wf','names','onOff','dates')
clear all

figure
for j=1:8
    subplot(2,4,j)
    for i=(j-1)*9+1:j*9
        tmp=wf(:,onOff<0,i);
        tmp=tmp-repmat(mean(tmp(50:450,:)),4500,1);
        plot(mean(tmp,2))
        hold on
    end
end