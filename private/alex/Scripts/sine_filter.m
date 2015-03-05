cd('/mnt/muench_data/user/alexandra/scripts')
%% Full Field Flicker Filter SINE
clear
cd('/mnt/muench_data/user/alexandra/scripts')
date='20130301_2';
words={'HL10','H30s','L30s'};
for wordsCNT=words
    codeWord=wordsCNT{1}
    
    
    mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/'];
    hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
    heka=dir([hekapath,'*.phys']);
    
    if ~exist([mainpath 'LinearFilters/'],'dir')
        mkdir([mainpath, 'LinearFilters']);
    end
    path2save=[mainpath 'LinearFilters/'];
    
    % make list of files to put into easy formatted unit (choosing by the name of the stimulus in heka files)
    file_list=[];
    for i=1:length(heka)
        if ~isempty(regexp(heka(i).name,codeWord, 'once'))
            file_list=[file_list i];
        end
    end
    
    timingPath=[mainpath 'timing/'];
    timings=dir([timingPath,'*.mat']);
    
    %get list of units to be reformatted
    basic_format_units_list=dir([mainpath, 'units/*.mat']);
    
    % get protocol
    protocol=read_header_field_heka(hekapath, heka(file_list(1)).name, 'Stimulus Protocol');
    protocol(protocol(:,3)==0&protocol(:,4)==0,4)=-1;
    protocol=round(protocol(2:end,[1,4]));
    protocol_accum=zeros(size(protocol,1),2,length(file_list));
    abused_accum=zeros(size(protocol,1),length(file_list));
    for i=1:length(file_list)
        protocol=read_header_field_heka(hekapath, heka(file_list(i)).name, 'Stimulus Protocol');
        protocol(protocol(:,3)==0&protocol(:,4)==0,4)=-1;
        protocol=round(protocol(2:end,[1,4]));
        protocol_accum(:,:,i)=protocol;
        load([timingPath, timings(file_list(i)).name]);
        abused=abused(2:end,5)-abused(2,5);
        abused_accum(:,i)=abused(3:end)*1000;
    end
    
    save([path2save,date,'_',codeWord,'_protocols_raw'],'abused_accum','protocol_accum')
    
    % load([path2save,date,'_',codeWord,'_protocols_raw'])
    
    tmp=(0:size(abused_accum,1))'*(7.85/3600);
    tmp=repmat(tmp(1:end-1),1,size(abused_accum,2));
    abused_accum=round(abused_accum+tmp);
    
    flicker=zeros(abused_accum(end,1),length(file_list));
    startingPoints=zeros(1,length(file_list));
    maxLength=zeros(1,length(file_list));
    
    for i=1:length(file_list)
        flips=abused_accum(1:end-2,i);
        pxls=(protocol_accum(1:end-2,2,i)-30)/30;
        startingPoints(i)=flips(1);
        flips=flips-flips(1)+1;
        
        tmp=zeros(flips(end),1);
        pxls(pxls==0)=100;
        tmp(flips(2:end),1)=pxls(2:end);
        % fill stimulus down to ms with brightness values where spikes
        % happen
        droppedMax=max(diff(flips));
        for j=1:droppedMax+1
            subFrame=flips(2:end)-j;
            subFrame(subFrame<1)=[];
            subFrame(find(tmp(subFrame,1)))=[];
            tmp(subFrame,1)=tmp(subFrame+1,1);
        end
        tmp(tmp(:,1)==100,1)=0;
        flicker(1:length(tmp),i)=tmp;
        maxLength(i)=length(tmp);
    end
    flicker=flicker(1:max(maxLength),:);
    
    save([path2save,date,'_flicker_',codeWord],'flicker','startingPoints','maxLength','abused_accum')
end

%% Full filter
filter_length=500;

names=cell(1,length(basic_format_units_list));
LinearFilter=zeros(filter_length,length(file_list),length(basic_format_units_list));
SpikeCount=zeros(length(file_list),length(basic_format_units_list));
FiringRate=zeros(abused_accum(end)+1000,length(file_list),length(basic_format_units_list));

for cnt=1:length(basic_format_units_list)
    basic_format_units_list(cnt).name
    load([mainpath,'units/',basic_format_units_list(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            a=convolved(spikes,40,0);
            FiringRate(1:length(a)-240,i,cnt)=a(121:end-120);
            spikes=spikes-startingPoints(i)+1;
            spikes(spikes<filter_length|spikes>size(flicker,1))=[];
            n=zeros(length(spikes),filter_length);
            for k=1:filter_length
                n(:,k)=flicker(spikes-k+1,i);
            end
            LinearFilter(1:filter_length,i,cnt)=sum(n);
        end
            SpikeCount(i,cnt)=length(spikes);
            names{cnt}=basic_format_units_list(cnt).name;

    end
end
save([path2save,date,'_LF_',codeWord],'LinearFilter','SpikeCount','names')
save([path2save,date,'_LF_',codeWord,'_FiringRate'],'FiringRate','names')

%% Plot Firing Rate
acc_tmp=0;
aver=2000;
for i=1:size(LinearFilter,4)
    tmp=FiringRate(:,3,i);
    t=floor(length(tmp)/aver);
    tmp=reshape(tmp(1:t*aver),aver,t);
    tmp=mean(tmp);
    acc_tmp=acc_tmp+tmp;
    plot(tmp)
    hold on
    for k=1:dim
        tt=round(abused_accum(per*(k-1)+1,3)/aver);
        line([tt,tt],[0,50],'color','r')
    end

end
plot(acc_tmp/size(LinearFilter,4),'color','k','linewidth',3)
k=1
round(abused_accum(per*(k-1)+1,3)/1000)

%% Plot

clear
date='20130301_1';
nd='87654321';
codeWord='sine';
load(['/mnt/muench_data/data/alexandra/MEA_data/',date,'/LinearFilters/',date,'_LF_',codeWord]);

path2save=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/',codeWord,'/'];

if ~exist(path2save,'dir')    
    mkdir(path2save);
end


figure(1)
set(gcf,'Position', [1 352 1011 625])
[rows,cols]=opt_subplots(length(nd));

for i=1:size(LinearFilter,3)
    tmp=LinearFilter(:,:,i);
    minn=min(tmp(:));
    maxx=max(tmp(:));
    for j=1:size(tmp,2)
        subplot(rows,cols,j);
        plot(tmp(:,j))
        line([0 500],[0 0],'color','k')
        title(['ND',nd(j)])
        axis([0 500 minn maxx])
    end
    subplot('Position',[0.5 0.97 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([names{i}(1:end-4),'_',codeWord],'FontSize',12,'FontWeight','bold','Interpreter','None')

    saveas(gcf,[path2save,names{i}(1:end-4),'_',codeWord,'.png'])
end



%% Partial filters: period
clear
date='20130301';
codeWord='sine';

per=3770; % period in frames

per=1885/5; % period in frames

mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/'];
basic_format_units_list=dir([mainpath, 'units/*.mat']);
path2take=[mainpath,'/LinearFilters/'];
load([path2take,date,'_flicker_',codeWord])
dim=floor((size(abused_accum,1)-2)/per);
filter_length=500;

hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);
file_list=[];
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,codeWord, 'once'))
        file_list=[file_list i];
    end
end

path2save=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/',codeWord,'_',int2str(per),'/'];

if ~exist(path2save,'dir')    
    mkdir(path2save);
end

names=cell(1,length(basic_format_units_list));
LinearFilter=zeros(filter_length,dim,length(file_list),length(basic_format_units_list));
SpikeCount=zeros(length(file_list),length(basic_format_units_list),dim);

for cnt=1:length(basic_format_units_list)
    basic_format_units_list(cnt).name
    load([mainpath,'units/',basic_format_units_list(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            spikes=spikes-startingPoints(i)+1;
            spikes(spikes<filter_length|spikes>size(flicker,1))=[];
            for t=1:dim
                spikesTMP=spikes(spikes>(abused_accum(per*(t-1)+1,i)-startingPoints(i))&spikes<(abused_accum(per*t,i)-startingPoints(i)));
                n=zeros(length(spikesTMP),filter_length);
                for k=1:filter_length
                    n(:,k)=flicker(spikesTMP-k+1,i);
                end
                LinearFilter(1:filter_length,t,i,cnt)=sum(n);
                SpikeCount(i,cnt,t)=length(spikesTMP);
            end
        end
        names{cnt}=basic_format_units_list(cnt).name;
    end
end
save([path2save,date,'_LF_',codeWord,'_',int2str(per)],'LinearFilter','SpikeCount','names')



 
 
 %% Plot partial filters
% 
% clear
% date='20130301';
% codeWord='sine';
% load(['/mnt/muench_data/data/alexandra/MEA_data/',date,'/LinearFilters/',date,'_LF_',codeWord]);

path2save=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/',codeWord,'_',int2str(per),'/'];

if ~exist(path2save,'dir')    
    mkdir(path2save);
end
nd='87654321';

figure(3)
set(gcf,'Position', [1 38 1680 939])
[rows,cols]=opt_subplots(length(nd));

for i=1:size(LinearFilter,4)
    tmp=LinearFilter(:,:,:,i);
    minn=min(tmp(:));
    maxx=max(tmp(:));
    for j=1:size(tmp,3)
        subplot(rows,cols,j);
        plot(mean(tmp(:,[1 4 5 8 9 12 13 16 17],j),2))
        hold on
        plot(mean(tmp(:,[2 3 6 7 10 11 14 15 18 19],j),2),'r')
        line([0 500],[0 0],'color','k')
        title(['ND',nd(j)])
        axis tight
        hold off
    end
    subplot('Position',[0.5 0.97 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([names{i}(1:end-4),'_',codeWord,'  Period ',int2str(per),' frames'],'FontSize',12,'FontWeight','bold','Interpreter','None')

    saveas(gcf,[path2save,names{i}(1:end-4),'_',codeWord,'_',int2str(per),'.png'])
end




refPointIndex=[150 130 100 85 85 95 90];
figure
nd='87654321';
for m=1:size(LinearFilter,4)
    for k=2:8
        subplot(2,4,k-1)
        tmp=LinearFilter(:,:,k,m);
%         if k<5
%             startInd=120;
%         else
%             startInd=90;
%         end
        for j=1:19
            b=mean(tmp(:,j:20:end),2);
%             sig(j,k)=mean(sigALL(j:20:end));
            plot(b,'color',[1 0 0]*j/19)
%             ind(j,k)=find(b(startInd:end)>0,1)+startInd;
            hold on
        end
        line([0 500], [0 0],'color','k')
        Ylimit=get(gca,'YLim');
        line([refPointIndex(k-1) refPointIndex(k-1)],[Ylimit(1) Ylimit(2)],'color','k')
        hold off
        title(['ND',nd(k)])
    end
    subplot('Position',[0.5 0.97 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([names{m}(1:end-4),'_',codeWord,'  Period ',int2str(per),' frames'],'FontSize',12,'FontWeight','bold','Interpreter','None')

    saveas(gcf,[path2save,names{m}(1:end-4),'_',codeWord,'_',int2str(per),'.png'])
end


onOff=zeros(1,size(LinearFilter,4))-1;
minimax=onOff;
for m=1:size(LinearFilter,4)    
    tmp=LinearFilter(80:180,:,3,m);
    tmp=mean(tmp(:,10:20:end),2);
    [~, b]=min(tmp);
    [~, c]=max(tmp);
    if b>c % on cell
        onOff(m)=1;
        minimax(m)=c+79;
    else
        minimax(m)=b+79;
    end
end

clear ind
for m=1:size(LinearFilter,4)    
    for k=2:8
        tmp=LinearFilter(:,:,k,m);
        
        for j=1:19
            b=mean(tmp(:,j:20:end),2);
            if onOff(m)>0
                [~, minPos]=max(b(1:150));
            else
                [~, minPos]=min(b(1:150));
            end
            temp=find(onOff(m)*b(minPos:end)<0,1)+minPos-1;
            if isempty(temp)
                ind(j,k-1,m)=0;
            else
                ind(j,k-1,m)=find(onOff(m)*b(minPos:end)<0,1)+minPos-1;
            end
        end
    end
end


clear sigALL
protocol=protocol_accum(:,2,1);
for i=1:dim
    sigALL(i)=std(protocol((i-1)*per+1:i*per));
end
clear sig
for j=1:19
    sig(j)=mean(sigALL(j:20:end));
end
sig=repmat(sig',1,7);

figure
for m=1:size(LinearFilter,4) 
    plot(ind(:,:,m),sig,'*');
    hold on
end


for k=1:7
    t=reshape(ind(:,k,:),19,size(LinearFilter,4) );
    for j=1:19
        meanInd(j,k)=mean(t(j,t(j,:)>70&t(j,:)<250));
    end
end
figure
plot(meanInd,sig,'*')

contr=sig(:,1)>3;
for m=1:size(LinearFilter,4)
    t=ind(contr,:,m);
    for k=1:7
        fitobject = fit(sig(contr,k),t(:,k),'exp1');
        p1(k,m)=fitobject.a;
        p2(k,m)=fitobject.b;
    end
end

p2(p2>0.02)=NaN;
p2(p2<-0.02)=NaN;

figure
plot(p2,'*')
hold on
plot(nanmean(p2'),'.','color','k','markersize',50)
set(gca,'xticklabel',{'ND7','ND6','ND5','ND4','ND3','ND2','ND1'})
line([1 7],[0,0],'color','k')



figure
a=nanmean(p1');
b=nanmean(p2');
for k=1:7
    plot(a(k).*sig(contr,k)+b(k),sig(contr,k),'color',[1 0 1]*k/7)
    hold on
end
legend({'ND7','ND6','ND5','ND4','ND3','ND2','ND1'})




p1(isnan(p2))=NaN;
figure
plot(p1)
hold on
plot(nanmean(p1'),'.','color','k','markersize',50)
set(gca,'xticklabel',{'ND7','ND6','ND5','ND4','ND3','ND2','ND1'})
line([1 7],[0,0],'color','k')





a=sig(:,6)
b=ind(:,6)

figure
j=3
tmp=LFtmp;
subplot(2,1,1)
 plot(tmp(:,:,j)) 
 hold on
 line([0 500],[0 0],'color','k')
 clear ind
 for i=1:size(LFtmp,2)
     ind(i)=find(tmp(120:end,i,j)>0,1);
 end
subplot(2,1,2)
 plot(ind,'*')
 b=diff(ind)
 mean(diff(ind))
 
 t=1;
 for i=1:10:length(ind)-5
     k(t)=ind(i+9)-ind(i);
     t=t+1;
 end
 
 


