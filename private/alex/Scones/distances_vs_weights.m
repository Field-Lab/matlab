%% datarun A

piece = '2013-10-10-5';
run = 'data002';
datarun = load_data(['/Volumes/Analysis/' piece '/streamed/' run '/' run]);
opt = struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',false);
datarun = load_data(datarun,opt);
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);
frame=6;
cellType=5;
load('/Volumes/Analysis/2013-10-10-5/stimuli/maps/map_data002_info.mat')

% plot STAs
a=zeros(length(info.cells), 1);
ncones=cell(length(info.cells), 1);
figure
for k=1:length(info.cells)
    a(k)=length(cones{k});
    if k>1
        tmp=cumsum(a);
        ncones{k}=(tmp(k-1)+1):tmp(k);
    else
        ncones{k}=1:a(k);
    end
    cellID=info.cells(k);
    cellInd=find(datarun.cell_ids==cellID);
    if isempty(datarun.stas.stas{cellInd})
        datarun = load_sta(datarun,'load_sta',cellID);
    end
    mySTA=squeeze(datarun.stas.stas{cellInd});
    mySTA=mySTA(:,:,frame);
    subplot(3,4,k)
    colormap gray
    imagesc(mySTA)
    hold on
    plot(stim.coord(ncones{k},1), stim.coord(ncones{k},2), 'rx')
    title(int2str(cellID))
    myWeights{k}=stim.weight(ncones{k});
end

% get distances
cone_distances=cell(length(info.cells),1);
nCentralCones=3;
allcones=false;
figure
for k=1:length(info.cells)
    cellID=info.cells(k);
    cellInd=find(datarun.cell_ids==cellID);
    mySTA=squeeze(datarun.stas.stas{cellInd});
    mySTA=mySTA(:,:,frame);
    if allcones
        mycones=1:length(stim.weight);
    else
        mycones=ncones{k};
    end
    coord=round(stim.coord(mycones, :));
    w=zeros(size(mycones));
    for i=1:length(mycones)
        w(i)=mySTA(coord(i,2),coord(i,1));
    end
    %         w=stim.weight(ncones{k});
    w=w/max(w);
    [w, ic]=sort(w, 'descend');
    
    x=stim.coord(mycones(ic),1)';
    y=stim.coord(mycones(ic),2)';
    
    x_center=x(1:nCentralCones);
    y_center=y(1:nCentralCones);
    w_center=w(1:nCentralCones);
    
    x_com=sum(w_center.*x_center)/sum(w_center);
    y_com=sum(w_center.*y_center)/sum(w_center);
    
    cone_distances{k}=sqrt((x-x_com).^2+(y-y_com).^2); % distances for current cell, ordered by weights
    subplot(3,4,k)
    plot(cone_distances{k}, w, '-x')
    
end

%% datarun B

piece = '2012-09-13-2';
run = 'data009';
datarun = load_data(['/Volumes/Analysis/' piece '/' run '/' run]);
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',false);
datarun = load_data(datarun,opt);
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);

cellType=5;
load('/Volumes/Analysis/2012-09-13-2/stimuli/maps/map_data009_info.mat')

% plot STAs
a=zeros(length(info.cells), 1);
ncones=cell(length(info.cells), 1);
figure
for k=1:length(info.cells)
    a(k)=length(cones{k});
    if k>1
        tmp=cumsum(a);
        ncones{k}=(tmp(k-1)+1):tmp(k);
    else
        ncones{k}=1:a(k);
    end
    cellID=info.cells(k);
    cellInd=find(datarun.cell_ids==cellID);
    if isempty(datarun.stas.stas{cellInd})
        datarun = load_sta(datarun,'load_sta',cellID);
    end
    mySTA=squeeze(datarun.stas.stas{cellInd});
    mySTA=mySTA(:,:,5);
    subplot(3,4,k)
    colormap gray
    imagesc(mySTA)
    hold on
    plot(stim.coord(ncones{k},1), stim.coord(ncones{k},2), 'rx')
    title(int2str(cellID))
    myWeights{k}=stim.weight(ncones{k});
end

% get distances
cone_distances=cell(length(info.cells),1);
nCentralCones=3;
allcones=true;
figure
for k=1:length(info.cells)
    cellID=info.cells(k);
    cellInd=find(datarun.cell_ids==cellID);
    mySTA=squeeze(datarun.stas.stas{cellInd});
    mySTA=mySTA(:,:,5);
    if allcones
        mycones=1:length(stim.weight);
    else
        mycones=ncones{k};
    end
    coord=round(stim.coord(mycones, :));
    w=zeros(size(mycones));
    for i=1:length(mycones)
        w(i)=mySTA(coord(i,2),coord(i,1));
    end
    %         w=stim.weight(ncones{k}); % weights for current cell, sorted in descending order
    w=w/max(w);
    subs=w(ncones{k});    
    [w, ic]=sort(w, 'descend');
    
    x=stim.coord(mycones(ic),1)';
    y=stim.coord(mycones(ic),2)';
    
    x_center=x(1:nCentralCones);
    y_center=y(1:nCentralCones);
    w_center=w(1:nCentralCones);
    
    x_com=sum(w_center.*x_center)/sum(w_center);
    y_com=sum(w_center.*y_center)/sum(w_center);
    
    cone_distances{k}=sqrt((x-x_com).^2+(y-y_com).^2); % distances for current cell, ordered by weights
    subplot(3,4,k)
    plot(cone_distances{k}, w, '-x')
    if allcones
        hold on
        tmp=find(ismember(w,subs));
        plot(cone_distances{k}(tmp), w(tmp), '-+r')
    end
    
end

% plot single cones for 1 cell

for k=1:length(info.cells)
    figure
    cellID=info.cells(k);
    cellInd=find(datarun.cell_ids==cellID);
    mySTA=squeeze(datarun.stas.stas{cellInd});
    mySTA=mySTA(:,:,5);
    
    coord=round(stim.coord);
    w=zeros(1,size(coord,1));
    for i=1:length(coord)
        w(i)=mySTA(coord(i,2),coord(i,1));
    end
    
    thr=sort(w);
    thr=thr(end-16);
    cnt=1;
    for i=1:length(coord)
        if mySTA(coord(i,2),coord(i,1))>thr
            subplot(4,4,cnt)
            colormap gray
            imagesc(mySTA(coord(i,2)-5:coord(i,2)+5, coord(i,1)-5:coord(i,1)+5))
            cnt=cnt+1;
        end
    end

    %         w=stim.weight(ncones{k}); % weights for current cell, sorted in descending order
    w=w/max(w);
    subs=w(ncones{k});    
    [w, ic]=sort(w, 'descend');
    
    x=stim.coord(mycones(ic),1)';
    y=stim.coord(mycones(ic),2)';
    
    x_center=x(1:nCentralCones);
    y_center=y(1:nCentralCones);
    w_center=w(1:nCentralCones);
    
    x_com=sum(w_center.*x_center)/sum(w_center);
    y_com=sum(w_center.*y_center)/sum(w_center);
    
    cone_distances{k}=sqrt((x-x_com).^2+(y-y_com).^2); % distances for current cell, ordered by weights
    subplot(3,4,k)
    plot(cone_distances{k}, w, '-x')
    if allcones
        hold on
        tmp=find(ismember(w,subs));
        plot(cone_distances{k}(tmp), w(tmp), '-+r')
    end
    
end

%% datarun D


piece = '2010-03-05-2';
run = 'data001';

% define data path
datarun = load_data(['/Volumes/Analysis/' piece '/' run '/' run]);
opt = struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);
datarun = load_data(datarun,opt);
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);

frame=5;
cellType=5;
load('/Volumes/Analysis/2010-03-05-2/stimuli/maps/map_data001_info.mat')

% plot STAs
a=zeros(length(info.cells), 1);
ncones=cell(length(info.cells), 1);
figure
for k=1:length(info.cells)
    a(k)=length(cones{k});
    if k>1
        tmp=cumsum(a);
        ncones{k}=(tmp(k-1)+1):tmp(k);
    else
        ncones{k}=1:a(k);
    end
    cellID=info.cells(k);
    cellInd=find(datarun.cell_ids==cellID);
    if isempty(datarun.stas.stas{cellInd})
        datarun = load_sta(datarun,'load_sta',cellID);
    end
    mySTA=squeeze(datarun.stas.stas{cellInd});
    mySTA=mySTA(:,:,frame);
    subplot(3,4,k)
    colormap gray
    imagesc(mySTA)
    hold on
     plot(stim.coord(ncones{k},1), stim.coord(ncones{k},2), 'rx')
    title(int2str(cellID))
    myWeights{k}=stim.weight(ncones{k});
end

% get distances
cone_distances=cell(length(info.cells),1);
nCentralCones=3;
allcones=true;
figure
for k=1:length(info.cells)
    cellID=info.cells(k);
    cellInd=find(datarun.cell_ids==cellID);
    mySTA=squeeze(datarun.stas.stas{cellInd});
    mySTA=mySTA(:,:,frame);
    if allcones
        mycones=1:length(stim.weight);
    else
        mycones=ncones{k};
    end
    coord=round(stim.coord(mycones, :));
    w=zeros(size(mycones));
    for i=1:length(mycones)
        w(i)=mySTA(coord(i,2),coord(i,1));
    end
    %         w=stim.weight(ncones{k});
    w=w/max(w);
    [w, ic]=sort(w, 'descend');
    
    x=stim.coord(mycones(ic),1)';
    y=stim.coord(mycones(ic),2)';
    
    x_center=x(1:nCentralCones);
    y_center=y(1:nCentralCones);
    w_center=w(1:nCentralCones);
    
    x_com=sum(w_center.*x_center)/sum(w_center);
    y_com=sum(w_center.*y_center)/sum(w_center);
    
    cone_distances{k}=sqrt((x-x_com).^2+(y-y_com).^2); % distances for current cell, ordered by weights
    subplot(3,4,k)
    plot(cone_distances{k}, w, '-x')
    
end

%% datarun E

piece = '2013-08-19-4';
run = 'data001';

% define data path
datarun = load_data(['/Volumes/Analysis/' piece '/streamed/' run '/' run]);
opt = struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);
datarun = load_data(datarun,opt);
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);


frame=5;
cellType=5;
load('/Volumes/Analysis/2013-08-19-4/stimuli/maps/map_data001_info.mat')

% plot STAs
a=zeros(length(info.cells), 1);
ncones=cell(length(info.cells), 1);
figure
for k=1:length(info.cells)
    a(k)=length(cones{k});
    if k>1
        tmp=cumsum(a);
        ncones{k}=(tmp(k-1)+1):tmp(k);
    else
        ncones{k}=1:a(k);
    end
    cellID=info.cells(k);
    cellInd=find(datarun.cell_ids==cellID);
    if isempty(datarun.stas.stas{cellInd})
        datarun = load_sta(datarun,'load_sta',cellID);
    end
    mySTA=squeeze(datarun.stas.stas{cellInd});
    mySTA=mySTA(:,:,frame);
    subplot(3,4,k)
    colormap gray
    imagesc(mySTA)
    hold on
    plot(stim.coord(ncones{k},1), stim.coord(ncones{k},2), 'rx')
    title(int2str(cellID))
    myWeights{k}=stim.weight(ncones{k});
end

% get distances
cone_distances=cell(length(info.cells),1);
nCentralCones=3;
allcones=true;
figure
for k=1:length(info.cells)
    cellID=info.cells(k);
    cellInd=find(datarun.cell_ids==cellID);
    mySTA=squeeze(datarun.stas.stas{cellInd});
    mySTA=mySTA(:,:,frame);
    if allcones
        mycones=1:length(stim.weight);
    else
        mycones=ncones{k};
    end
    coord=round(stim.coord(mycones, :));
    w=zeros(size(mycones));
    for i=1:length(mycones)
        w(i)=mySTA(coord(i,2),coord(i,1));
    end
    %         w=stim.weight(ncones{k});
    w=w/max(w);
    [w, ic]=sort(w, 'descend');
    
    x=stim.coord(mycones(ic),1)';
    y=stim.coord(mycones(ic),2)';
    
    x_center=x(1:nCentralCones);
    y_center=y(1:nCentralCones);
    w_center=w(1:nCentralCones);
    
    x_com=sum(w_center.*x_center)/sum(w_center);
    y_com=sum(w_center.*y_center)/sum(w_center);
    
    cone_distances{k}=sqrt((x-x_com).^2+(y-y_com).^2); % distances for current cell, ordered by weights
    subplot(3,4,k)
    plot(cone_distances{k}, w, '-x')
    
end

%% datarun H

piece = '2013-08-19-5';
run = 'data001';

% define data path
datarun = load_data(['/Volumes/Analysis/' piece '/streamed/' run '/' run]);
opt = struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);
datarun = load_data(datarun,opt);
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);



frame=5;
cellType=5;
load('/Volumes/Analysis/2013-08-19-5/stimuli/maps/map_data001_info.mat')

% plot STAs
a=zeros(length(info.cells), 1);
ncones=cell(length(info.cells), 1);
figure
for k=1:length(info.cells)
    a(k)=length(cones{k});
    if k>1
        tmp=cumsum(a);
        ncones{k}=(tmp(k-1)+1):tmp(k);
    else
        ncones{k}=1:a(k);
    end
    cellID=info.cells(k);
    cellInd=find(datarun.cell_ids==cellID);
    if isempty(datarun.stas.stas{cellInd})
        datarun = load_sta(datarun,'load_sta',cellID);
    end
    mySTA=squeeze(datarun.stas.stas{cellInd});
    mySTA=mySTA(:,:,frame);
    subplot(3,4,k)
    colormap gray
    imagesc(mySTA)
    hold on
    plot(stim.coord(ncones{k},1), stim.coord(ncones{k},2), 'rx')
    title(int2str(cellID))
    myWeights{k}=stim.weight(ncones{k});
end

% get distances
cone_distances=cell(length(info.cells),1);
nCentralCones=3;
allcones=true;
figure
for k=1:length(info.cells)
    cellID=info.cells(k);
    cellInd=find(datarun.cell_ids==cellID);
    mySTA=squeeze(datarun.stas.stas{cellInd});
    mySTA=mySTA(:,:,frame);
    if allcones
        mycones=1:length(stim.weight);
    else
        mycones=ncones{k};
    end
    coord=round(stim.coord(mycones, :));
    w=zeros(size(mycones));
    for i=1:length(mycones)
        w(i)=mySTA(coord(i,2),coord(i,1));
    end
    %         w=stim.weight(ncones{k});
    w=w/max(w);
    [w, ic]=sort(w, 'descend');
    
    x=stim.coord(mycones(ic),1)';
    y=stim.coord(mycones(ic),2)';
    
    x_center=x(1:nCentralCones);
    y_center=y(1:nCentralCones);
    w_center=w(1:nCentralCones);
    
    x_com=sum(w_center.*x_center)/sum(w_center);
    y_com=sum(w_center.*y_center)/sum(w_center);
    
    cone_distances{k}=sqrt((x-x_com).^2+(y-y_com).^2); % distances for current cell, ordered by weights
    subplot(3,4,k)
    plot(cone_distances{k}, w, '-x')
    
end
