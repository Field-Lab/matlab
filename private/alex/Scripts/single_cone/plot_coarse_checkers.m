% load a CBF filter. Run from windows
path2data='S:\user\alexandra\MEA_data\20121024\CBF_filter\';
units=dir([path2data,'*CBFlicker_LF*.mat']);
load([path2data,units(1).name])

% get HeatMaps
% close all
m=figure(3);
set(m,'Position',[1 31 1680 946])
pos=[0.01,0.75, 0.93  0.12];
timePoint2show=[20:20:500];
cnt=0;
all_acc=0;
for trial=1:3   
    tmp=CBFfilter(:,:,trial);
    tmp=tmp/max(abs(tmp(:)));
    b=mean(tmp);
    acc=[];
    for i=1:length(timePoint2show)
        a=reshape(tmp(:,timePoint2show(i)),40,40);
        a=flipud(a);
        a=fliplr(a);
        acc=[acc a ones(40,1)];
    end
    acc(:,end-2:end)=[];
    all_acc=all_acc+acc;
    subplot(3,1,trial)
    colormap('gray')
    imagesc(all_acc)
%     MyHeatMap(all_acc/(cnt+1),m,k,1);
    cnt=cnt+1;
end
for j=1:length(timePoint2show)
    text(15+(j-1)*33,32,int2str(500-timePoint2show(j)),'Fontsize',8,'FontWeight','bold','horizontalAlignment','center')
end


% get and plot filter as a single line - full 8min sequences

close all
m=figure(4);
set(m,'Position',[1 31 1680 946])

path2data='S:\user\alexandra\MEA_data\20121024\CBF_filter\';
units=dir([path2data,'*CBFlicker_LF*.mat']);
load([path2data,units(1).name])

findHotSpot=1; % trial number
for units_cnt=1:length(units)
    load([path2data,units(units_cnt).name])    
  
    tmp=CBFfilter(:,:,findHotSpot);
    tmp=tmp';
    [val co]=max(std(tmp));
    centrSize=2;
    surrSize=7;
    b=tmp(:,co);
    centr=[co-centrSize:co+centrSize];
    i=[];
    for j=-centrSize:centrSize
        i=[i centr+40*j];
    end
    i(i>1600)=[];
    i(i<1)=[];
    
    
    acc=zeros(1,500);
    for trial=1:3
        tmp=CBFfilter(:,:,trial);
        tmp=tmp/max(abs(tmp(:)));
        acc(trial,1:500)=mean(tmp(i,:));
%         acc(trial,1:500)=acc(trial,1:500)/max(abs(acc(trial,1:500)));
        acc(trial,1:500)=acc(trial,500:-1:1);

    end
    a=regexp(units(units_cnt).name,'unit');
    tit=units(units_cnt).name(a+5:a+8);
    subplot(6,7,units_cnt)
    acc=acc(:,1:500)';
    plot(acc)
    title(tit)
    axis tight
end




load('S:\user\alexandra\MEA_data\20121024\cones\A20121024_CH10_sort1_100_unit_0008_cones.mat')
m=figure(3);
set(m,'Position',[1 31 1680 946])
pos=[0.01,0.75, 0.93  0.12];
timePoint2show=[20:20:500];
cnt=0;
all_acc=0;
tmp=coneFilter;
coneFilter=coneFilter/max(abs(coneFilter(:)));
b=mean(coneFilter);
acc=[];
for i=1:length(timePoint2show)
    a=reshape(coneFilter(:,timePoint2show(i)),50,50);
    acc=[acc a ones(50,1)];
end
acc(:,end-2:end)=[];
all_acc=all_acc+acc;
colormap('gray')
imagesc(all_acc)