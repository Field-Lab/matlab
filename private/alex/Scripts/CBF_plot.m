% load a CBF filter. Run from windows
path2data='S:\user\alexandra\MEA_data\20121024\easy_formatted_units\';
units=dir([path2data,'*spike_info_CBFlicker*.mat']);
load([path2data,units(1).name])

% get HeatMaps
close all
m=figure(3);
set(m,'Position',[1 31 1680 946])
pos=[0.01,0.75, 0.93  0.12];
timePoint2show=[700:20:950];
cnt=0;
all_acc=0;
for trial=3:9    
    tmp=CBFfilter(:,:,trial);
    tmp=tmp/max(abs(tmp(:)));
    b=mean(tmp);
    acc=[];
    for i=1:length(timePoint2show)
        a=reshape(tmp(:,timePoint2show(i)),40,40);
        a=flipud(a);
        a=fliplr(a);
        acc=[acc a ones(40,3)];
    end
    acc(:,end-2:end)=[];
    all_acc=all_acc+acc;
    k=subplot('Position',[pos(1) pos(2)-0.12*cnt pos(3:4)]);
    MyHeatMap(all_acc/(cnt+1),m,k,1);
    cnt=cnt+1;
end
for j=1:length(timePoint2show)
    text(15+(j-1)*33,32,int2str(1000-timePoint2show(j)),'Fontsize',8,'FontWeight','bold','horizontalAlignment','center')
end


% get and plot filter as a single line - full 8min sequences

close all
m=figure(1);
set(m,'Position',[1 31 1680 946])

filterspath='S:\user\alexandra\MEA_data\20120628\CBF_1000\';
units=dir([filterspath,'*.mat']);
findHotSpot=1; % trial number
for units_cnt=1:length(units)
    load([filterspath,units(units_cnt).name])    
  
    tmp=CBFfilter(:,:,findHotSpot);
    tmp=tmp';
    [val co]=max(std(tmp));
    centrSize=3;
    surrSize=7;
    b=tmp(:,co);
    centr=[co-centrSize:co+centrSize];
    i=[];
    for j=-centrSize:centrSize
        i=[i centr+40*j];
    end
    i(i>1600)=[];
    i(i<1)=[];
    
    
    acc=zeros(1,1000);
    for trial=1:3
        tmp=CBFfilter(:,:,trial);
        tmp=tmp/max(abs(tmp(:)));
        acc(trial,1:1000)=mean(tmp(i,:));
        acc(trial,1:1000)=acc(trial,1:1000)/max(abs(acc(trial,1:1000)));
        acc(trial,1:1000)=acc(trial,1000:-1:1);

    end
    a=regexp(units(units_cnt).name,'unit');
    tit=units(units_cnt).name(a+5:a+8);
    subplot(4,6,units_cnt)
    acc=acc(:,1:500)';
    plot(acc)
    title(tit)
    axis tight
end

figure(1)
for i=1:3
    subplot(2,3,i)
    a=std(CBFfilter(:,:,i),0,2);
    imagesc(reshape(a,40,40));
    title(['ND',int2str(5-i)])
    subplot(2,3,i+3)
    plot(acc(:,i),'linewidth',2)
    title(['ND',int2str(5-i)])
    line([100,100],[-1,1],'color','k')
end
% get and plot filter as a single line - short sequences

close all
m=figure(4);
set(m,'Position',[1 31 1680 946])
filterspath='S:\user\alexandra\MEA_data\20120628\CBF_1000_short\';
units=dir([filterspath,'*60s.mat']);
for units_cnt=1:length(units)
    load([filterspath,units(units_cnt).name])    
  
    tmp=CBFfilter(:,:,3);
    tmp=tmp';
    [val co]=max(std(tmp));
    centrSize=3;
    surrSize=7;
    b=tmp(:,co);
    centr=[co-centrSize:co+centrSize];
    i=[];
    for j=-centrSize:centrSize
        i=[i centr+40*j];
    end
    i(i>1600)=[];
    i(i<1)=[];
    
    
    acc=zeros(1,1000);
    for trial=1:9
        tmp=CBFfilter(:,:,trial);
        tmp=tmp/max(abs(tmp(:)));
        acc(trial,1:1000)=mean(tmp(i,:));
        acc(trial,1:1000)=acc(trial,1:1000)/max(abs(acc(trial,1:1000)));
        acc(trial,1:1000)=acc(trial,1000:-1:1);

    end
    a=regexp(units(units_cnt).name,'unit');
    tit=units(units_cnt).name(a+5:a+8);
    subplot(4,6,units_cnt)
    acc=acc(:,1:500)';
    plot(acc)
    title(tit)
    axis tight
end
