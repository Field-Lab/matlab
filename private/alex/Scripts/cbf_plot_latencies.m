clear
path2filter=['S:\data\alexandra\MEA_data\20120902_2\CBF_filter\'];
filters=dir([path2filter,'\*.mat']);
findHotSpot=1; % trial number
nds='432'
for unit=1:length(filters)
    load([path2filter,filters(unit).name]);
    name=filters(unit).name(1:end-17);    
    
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
    
    acc=zeros(3,500);
    for trial=1:3
        tmp=CBFfilter(:,:,trial);
        tmp=tmp/max(abs(tmp(:)));
        acc(trial,1:500)=mean(tmp(i,:));
        acc(trial,1:500)=acc(trial,1:500)/max(abs(acc(trial,1:500)));
        acc(trial,1:500)=acc(trial,500:-1:1);        
    end
    acc=acc';
    figure
    cnt=1;
    for i=1:3
        subplot(1,3,cnt)
        plot(acc(:,i))
        [~,k(1)]=min(acc(:,i));
        [~,k(2)]=max(acc(:,i));
        title(['ND',nds(cnt),', ',int2str(SpikeCount(i)),'sp, ',int2str(k)])
        axis tight
        cnt=cnt+1;
    end
    saveas(gcf,['S:\data\alexandra\MEA_data\20120902_2\CBFplot\',name,'_CBFplot.emf'])
    close all
end



clear
path2filter=['S:\data\alexandra\MEA_data\20120902_1\CBF_filter\'];
filters=dir([path2filter,'\*.mat']);
findHotSpot=1; % trial number
nds='432'
figure

for unit=1:length(filters)
    load([path2filter,filters(unit).name]);
    name=filters(unit).name(1:end-17);    
    
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
    
    acc=zeros(3,500);
    for trial=1:3
        tmp=CBFfilter(:,:,trial);
        tmp=tmp/max(abs(tmp(:)));
        acc(trial,1:500)=mean(tmp(i,:));
        acc(trial,1:500)=acc(trial,1:500)/max(abs(acc(trial,1:500)));
        acc(trial,1:500)=acc(trial,500:-1:1);        
    end
    acc=acc';
    
    subplot(5,6,unit)
    plot(acc(1:250,1:3))
    title(name,'interpreter','none')
    axis tight
    pause(3)
end
% saveas(gcf,['S:\data\alexandra\MEA_data\20120902_1\CBFplot\',name,'_CBFplot.emf'])
% close all