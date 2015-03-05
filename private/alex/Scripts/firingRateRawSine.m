
clear
dates=cell(1,3);
dates{1}='20130301';
dates{2}='20130301_1';
dates{3}='20130301_2';

codeWord='HL10';
per=600; % period in frames (start period). For sine:130, for HL10: 600, for H30s and L30s: 1800
dur=per*1; % duration of interval. for sine: per*5, for others: per*1
cntDate=1;
names=cell(1,103);
FiringRate=zeros(11000,16,103,6);
isi=zeros(100,16,103,6);
path2save=['/mnt/muench_data/data/alexandra/MEA_data/sine_analysis/HL_raw_firingRate/'];

if ~exist(path2save,'dir')
    mkdir(path2save);
end
for ttt=dates
    date=cell2mat(ttt)
    mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/'];
    basic_format_units_list=dir([mainpath, 'units/*.mat']);
    path2take=[mainpath,'LinearFilters/'];
    load([path2take,date,'_flicker_',codeWord])
    dim=floor((size(abused_accum,1)-2)/per);
    while (size(abused_accum,1)-2)<(per*(dim-1)+dur)
        dim=dim-1;
    end
    
    hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
    heka=dir([hekapath,'*.phys']);
    file_list=[];
    for i=1:length(heka)
        if ~isempty(regexp(heka(i).name,codeWord, 'once'))
            file_list=[file_list i];
        end
    end
    
    
    for cnt=1:length(basic_format_units_list)

        load([mainpath,'units/',basic_format_units_list(cnt).name]);
        for i=1:length(file_list)
            spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
            if ~isempty(spikes)
                spikes=spikes-startingPoints(i)+1;                
                for t=1:dim
                    spikesTMP=spikes(spikes>(abused_accum(per*(t-1)+1,i)-startingPoints(i))&spikes<(abused_accum(per*(t-1)+dur,i)-startingPoints(i)));
                    if t>1
                        spikesTMP=spikesTMP-(abused_accum(per*(t-2)+dur,i)-startingPoints(i));
                    end
                    tmp=convolved(spikesTMP,40,11000);
                    FiringRate(:,i,cntDate,t)=tmp(121:end-120);
                    a=diff(spikesTMP);
                    a(a>100)=[];
                    if length(a)>4
                        a(a==0)=1;
                        [z,x]=hist(a,unique(a));
                        isi(x,i,cntDate,t)=isi(x,i,cntDate,t)+z';
                    end
                end
            end            
        end
        name=basic_format_units_list(cnt).name(1:end-4);
        if cntDate>=33&&cntDate<67
            name=[name(1:9),'_1',name(10:end)];
        elseif cntDate>=67
            name=[name(1:9),'_2',name(10:end)];
        end
        names{cntDate}=basic_format_units_list(cnt).name;
        cntDate=cntDate+1;
    end

end

save([path2save,'FiringRateRaw'],'isi','FiringRate','names')
   

%% Plots
codeWord='HL10';
load(['/mnt/muench_data/data/alexandra/MEA_data/sine_analysis/',codeWord,'_accumParameters'])
load(['/mnt/muench_data/data/alexandra/MEA_data/sine_analysis/HL_raw_firingRate/','FiringRateRaw'],'isi','FiringRate','names')
goodUnits=xlsread('/mnt/muench_data/data/alexandra/MEA_data/sine_analysis/ListOfGoodFiltersPerND_sine.xls');
onOff_103=goodUnits(:,9);
goodUnits=goodUnits(:,1:8);

figure
kk=1;
i=7;
for cnt=1:103
    if onOff_103(cnt)>0
        a=mean(FiringRate(:,i,cnt,1:2:end),4);
        subplot(5,5,kk)
        plot(a(1:10000))
        hold on
        b=mean(FiringRate(:,i,cnt,2:2:end),4);                
        plot(b(1:10000),'r')
        if mean(SpikeCount(i,cnt,2:2:end))>0
            c=(round(mean(SpikeCount(i,cnt,1:2:end))/mean(SpikeCount(i,cnt,2:2:end))*10))/10;
        end
        d=(round(corr(a(1:10000),b(1:10000))*100))/100;
        title([int2str(cnt),' high/low=',num2str(c),', corr=',num2str(d)])
        kk=kk+1;        
    end
end



%% GET FIRING RATE FOR SINE
clear
dates=cell(1,3);
dates{1}='20130301';
dates{2}='20130301_1';
dates{3}='20130301_2';

codeWord='sine';
cntDate=1;
names=cell(1,103);
FiringRate=zeros(610000,8,103);
path2save=['/mnt/muench_data/data/alexandra/MEA_data/sine_analysis/sine_raw_firingRate/'];

if ~exist(path2save,'dir')
    mkdir(path2save);
end
for ttt=dates
    date=cell2mat(ttt)
    mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/'];
    basic_format_units_list=dir([mainpath, 'units/*.mat']);
    path2take=[mainpath,'LinearFilters/'];
    load([path2take,date,'_flicker_',codeWord])
    
    hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
    heka=dir([hekapath,'*.phys']);
    file_list=[];
    for i=1:length(heka)
        if ~isempty(regexp(heka(i).name,codeWord, 'once'))
            file_list=[file_list i];
        end
    end
    
    
    for cnt=1:length(basic_format_units_list)

        load([mainpath,'units/',basic_format_units_list(cnt).name]);
        for i=1:length(file_list)
            spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
            if ~isempty(spikes)
                tmp=convolved(spikes,40,610000);
                FiringRate(:,i,cntDate)=tmp(121:end-120);
            end            
        end
        name=basic_format_units_list(cnt).name(1:end-4);
        if cntDate>=33&&cntDate<67
            name=[name(1:9),'_1',name(10:end)];
        elseif cntDate>=67
            name=[name(1:9),'_2',name(10:end)];
        end
        names{cntDate}=basic_format_units_list(cnt).name;
        cntDate=cntDate+1;
    end

end

save([path2save,'FiringRateRaw'],'FiringRate','names')
 
%% Plots
codeWord='sine';
load(['/mnt/muench_data/data/alexandra/MEA_data/sine_analysis/',codeWord,'_accumParameters'])
goodUnits=xlsread('/mnt/muench_data/data/alexandra/MEA_data/sine_analysis/ListOfGoodFiltersPerND_sine.xls');
onOff_103=goodUnits(:,9);
goodUnits=goodUnits(:,1:8);

figure
kk=1;
i=2;
for cnt=1:103
%     if onOff_103(cnt)>0
%         b=(round(correlFiringRate(i,cnt)*100))/100;
%         if b<-0.2
            subplot(8,6,kk)
            if onOff(cnt)>0
                plot(FiringRate(1:600000,i,cnt),'r')
            else
                plot(FiringRate(1:600000,i,cnt),'b')
            end
            b=(round(correlFiringRate(i,cnt)*100))/100;
            title([int2str(cnt),'  correlation=',num2str(b), ',   ',int2str(onOff(cnt))])
            kk=kk+1;
%         end
%     end
end

figure
plot(FiringRate(1:600000,i,22),'r')


i=3
cnt=15
a=zeros(2000,103);

figure
for i=1:8
    kk=1;
    for dim=1:300:600000
        a(kk,:)=reshape(std(FiringRate(dim:dim+300,i,:)),103,1);
        kk=kk+1;
    end
    subplot(4,2,i)
    hold on
    c=onOff>0&correlFiringRate(i,:)>0.6;
    plot(mean(a(:,c)'),'b')
    c=onOff>0&correlFiringRate(i,:)<-0.2;
    plot(mean(a(:,c)'),'r')
    title(['ND',int2str(9-i)])
    drawnow
end



b=corr(reshape(FiringRate(:,3,onOff<0),610000,sum(onOff<0)));
figure
imagesc(b)


b=corr(reshape(FiringRate(:,3,:),610000,103));
figure
imagesc(b)


%% GET FIRING RATE FOR HL10
clear
dates=cell(1,3);
dates{1}='20130301';
dates{2}='20130301_1';
dates{3}='20130301_2';

codeWord='HL10';
cntDate=1;
names=cell(1,103);
FiringRate=zeros(65000,16,103);
path2save=['/mnt/muench_data/data/alexandra/MEA_data/sine_analysis/sine_raw_firingRate/'];

if ~exist(path2save,'dir')
    mkdir(path2save);
end
for ttt=dates
    date=cell2mat(ttt)
    mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/'];
    basic_format_units_list=dir([mainpath, 'units/*.mat']);
    path2take=[mainpath,'LinearFilters/'];
    load([path2take,date,'_flicker_',codeWord])
    
    hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
    heka=dir([hekapath,'*.phys']);
    file_list=[];
    for i=1:length(heka)
        if ~isempty(regexp(heka(i).name,codeWord, 'once'))
            file_list=[file_list i];
        end
    end
    
    
    for cnt=1:length(basic_format_units_list)

        load([mainpath,'units/',basic_format_units_list(cnt).name]);
        for i=1:length(file_list)
            spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
            if ~isempty(spikes)
                tmp=convolved(spikes,40,65000);
                FiringRate(:,i,cntDate)=tmp(121:end-120);
            end            
        end
        name=basic_format_units_list(cnt).name(1:end-4);
        if cntDate>=33&&cntDate<67
            name=[name(1:9),'_1',name(10:end)];
        elseif cntDate>=67
            name=[name(1:9),'_2',name(10:end)];
        end
        names{cntDate}=basic_format_units_list(cnt).name;
        cntDate=cntDate+1;
    end

end

save([path2save,'FiringRateRaw_HL10'],'FiringRate','names')

 %% Plots
codeWord='HL10';
load(['/mnt/muench_data/data/alexandra/MEA_data/sine_analysis/',codeWord,'_accumParameters'])
goodUnits=xlsread('/mnt/muench_data/data/alexandra/MEA_data/sine_analysis/ListOfGoodFiltersPerND_sine.xls');
onOff_103=goodUnits(:,9);
goodUnits=goodUnits(:,1:8);

figure
kk=1;
i=3;
for cnt=31:103
%     if onOff_103(cnt)>0
%         b=(round(correlFiringRate(i,cnt)*100))/100;
%         if b<-0.2
            subplot(6,4,kk)
            if onOff(cnt)>0
                plot(FiringRate(:,i,cnt),'r')
            else
                plot(FiringRate(:,i,cnt),'b')
            end
            title([int2str(cnt),',   ',int2str(onOff(cnt))])
%             b=(round(correlFiringRate(i,cnt)*100))/100;
%             title([int2str(cnt),'  correlation=',num2str(b), ',   ',int2str(onOff(cnt))])
            kk=kk+1;
%         end
%     end
end
figure
  plot(FiringRate(:,i,48),'r')