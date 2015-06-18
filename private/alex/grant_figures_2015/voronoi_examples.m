% load datarun

datarun = load_data('/Volumes/Analysis/2010-03-05-2/d12-15-norefit/data015-from-d12_15/data015-from-d12_15');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);

vorrun = load_data('/Volumes/Analysis/2010-03-05-2/d12-15-norefit/data014-from-d12_15/data014-from-d12_15');
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_neurons(vorrun);
[inputs_v, refresh_v, duration_v] = get_wn_movie_ath(vorrun, 'BW-1-8-0.48-11111-1517x1.xml');
vormap = load('/Volumes/Archive/2010-03-05-2/Visual/cone_map-data-014.txt');

datarun1 = load_data('/Volumes/Analysis/2010-03-05-2/data013/data013');
datarun1 = load_params(datarun1,'verbose',1);
datarun1 = load_cones_ath(datarun1,'streamed');

% creat voronoi overview with rectangles for select regions

figsize = 135;
figure
set(gcf,'units', 'points','position',[400   400  figsize, figsize])
set(gca, 'visible', 'off', 'dataaspectratio',[1 1 1])
subplot('position',[0 0 1 1])
voronoi_plot_no_centers(datarun1.cones.centers(:,1), datarun1.cones.centers(:,2))
axis ij
hold on
axis([58 270 58 270])
rectangle('position',[60 90 116 116], 'Edgecolor','b','linewidth',2);
rectangle('position',[85 133 23 23], 'Edgecolor','r', 'linewidth',2);


%% voronoi sta params

sta_params.length = 6;
sta_params.offset = 0;

%% large cell

datarunID = 57;

my_sta=zeros(size(inputs_v,1),sta_params.length);

spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh_v)); % spikes in frames
spikes(spikes<sta_params.length-sta_params.offset)=[];

nspikes = length(spikes);

while ~isempty(spikes)
    [~, ia, ~] = unique(spikes);
    for j=1:sta_params.length
        my_sta(:,sta_params.length-j+1) = my_sta(:, sta_params.length-j+1) +...
            sum(inputs_v(:,spikes(ia)-sta_params.length+j+sta_params.offset),2);
    end
    spikes(ia)=[];
end
my_sta=my_sta/nspikes;

% get voronoi map sta
tmp_map = vormap;
tt=0:max(tmp_map(:));
vorsta=zeros(320,320,sta_params.length);
cnt = 1;
for i=1:1517
    [a, b] = find(tmp_map==tt(i+1));
    if ~isempty(a)
        for j = 1:length(a)
            vorsta(a(j),b(j),:) = my_sta(cnt,:);
        end
    end
    cnt=cnt+1;
end

% figure
% for i=1:6    
%     subplot(2,3,i)
%     colormap gray
%     imagesc(vorsta(:,:,i));
% end
% 
% figure
% vorsta_tmp = vorsta(:, :, 2);
% colormap gray
% imagesc(vorsta_tmp)

tmp = my_sta;
% figure
% plot(my_sta')

thr = robust_std(tmp(:,3))*2;
thr1 = robust_std(tmp(:,2))*2;

tc = mean(tmp(tmp(:,2)<-thr1 & tmp(:,3)>thr,:));

% tc = mean(tmp(tmp(:,2)<-0.0125 & tmp(:,3)>0.0125,:));
% figure
% plot(tc)
tmp = my_sta(:,1:6).*repmat(tc(:,1:6), size(my_sta,1),1);
tmp = sum(tmp,2);

% get voronoi map sta
tmp_map = vormap;
tt=0:max(tmp_map(:));
vorsta=zeros(320,320);
cnt = 1;
for i=1:1517
    [a, b] = find(tmp_map==tt(i+1));
    if ~isempty(a)
        for j = 1:length(a)
            vorsta(a(j),b(j)) = tmp(cnt);
        end
    end
    cnt=cnt+1;
end

%%%%%%%%%% PLOT LARGE CELL VORONOI %%%%%%%%%%
figsize = 135;
figure
set(gcf,'units', 'points','position',[400   400  figsize, figsize])
colormap gray
subplot('position',[0 0 1 1])
imagesc(-vorsta)
axis([60 176 90 206])
set(gca, 'visible', 'off', 'dataaspectratio',[1 1 1])

% STA 

sta = squeeze(double(datarun.stas.stas{datarunID}));

% figure
% for i=1:6    
%     subplot(2,3,i)
%     colormap gray
%     imagesc(sta(:,:,i));
% end
% 
% figure
% colormap gray
% imagesc(sta(:,:,4));

tmp = sta(:,:,4);
thr = robust_std(tmp(:))*2;
tmp1 = sta(:,:,3);
thr1 = robust_std(tmp1(:))*2;
coords =  find(tmp<-thr & tmp1>thr1);

tc = [];
for i=1:5
    tmp = sta(:,:,i);
    tc(i)=mean(tmp(coords));
end
% figure
% plot(tc)
tt = shiftdim(repmat(tc',1,320,320),1);
tmp = sta(:,:,1:5).*tt;
tmp = sum(tmp,3);


%%%%%%%%%% PLOT LARGE CELL STA %%%%%%%%%%
figsize = 135;
figure
set(gcf,'units', 'points','position',[400   400  figsize, figsize])
colormap gray
subplot('position',[0 0 1 1])
imagesc(-tmp)
axis([60 176 90 206])
set(gca, 'visible', 'off', 'dataaspectratio',[1 1 1])


%% ON parasol

datarunID = 107;

my_sta=zeros(size(inputs_v,1),sta_params.length);

spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh_v)); % spikes in frames
spikes(spikes<sta_params.length-sta_params.offset)=[];

nspikes = length(spikes);

while ~isempty(spikes)
    [~, ia, ~] = unique(spikes);
    for j=1:sta_params.length
        my_sta(:,sta_params.length-j+1) = my_sta(:, sta_params.length-j+1) +...
            sum(inputs_v(:,spikes(ia)-sta_params.length+j+sta_params.offset),2);
    end
    spikes(ia)=[];
end
my_sta=my_sta/nspikes;

% get voronoi map sta
tmp_map = vormap;
tt=0:max(tmp_map(:));
vorsta=zeros(320,320,sta_params.length);
cnt = 1;
for i=1:1517
    [a, b] = find(tmp_map==tt(i+1));
    if ~isempty(a)
        for j = 1:length(a)
            vorsta(a(j),b(j),:) = my_sta(cnt,:);
        end
    end
    cnt=cnt+1;
end

% figure
% for i=1:6    
%     subplot(2,3,i)
%     colormap gray
%     imagesc(vorsta(:,:,i));
% end
% 
% figure
% vorsta_tmp = vorsta(:, :, 2);
% colormap gray
% imagesc(vorsta_tmp)

tmp = my_sta;
% figure
% plot(my_sta')

thr = robust_std(tmp(:,3))*2;
thr1 = robust_std(tmp(:,2))*2;
tc = mean(tmp(tmp(:,2)>thr1 & tmp(:,3)<-thr,:)); % ON cells
% tc = mean(tmp(tmp(:,2)<-thr1 & tmp(:,3)>thr,:)); % OFF cells
% figure
% plot(tc)
tmp = my_sta(:,1:6).*repmat(tc(:,1:6), size(my_sta,1),1);
tmp = sum(tmp,2);

% get voronoi map sta
tmp_map = vormap;
tt=0:max(tmp_map(:));
vorsta=zeros(320,320);
vorsta_tmp = vorsta;
cnt = 1;
for i=1:1517
    [a, b] = find(tmp_map==tt(i+1));
    if ~isempty(a)
        for j = 1:length(a)
            vorsta(a(j),b(j)) = tmp(cnt);
            
            vorsta_tmp(a(j),b(j)) = my_sta(cnt,2);
        end
    end
    cnt=cnt+1;
end


%%%%%%%%%% PLOT PARASOL CELL VORONOI %%%%%%%%%%
figsize = 135;
figure
set(gcf,'units', 'points','position',[400   400  figsize, figsize])
colormap gray
subplot('position',[0 0 1 1])
imagesc(vorsta)
axis([60 176 90 206])
set(gca, 'visible', 'off', 'dataaspectratio',[1 1 1])


% STA 

sta = squeeze(double(datarun.stas.stas{datarunID}));
% 
% figure
% for i=1:6    
%     subplot(2,3,i)
%     colormap gray
%     imagesc(sta(:,:,i));
% end


tmp = sta(:,:,4);
thr = robust_std(tmp(:))*2;
tmp1 = sta(:,:,3);
thr1 = robust_std(tmp1(:))*2;

coords =  find(tmp>thr & tmp1<-thr1); % ON cells
% coords =  find(tmp<-thr & tmp1>thr1); % OFF cells

tc = [];
for i=1:5
    tmp = sta(:,:,i);
    tc(i)=mean(tmp(coords));
end
% figure
% plot(tc)
tt = shiftdim(repmat(tc',1,320,320),1);
tmp = sta(:,:,1:5).*tt;
tmp = sum(tmp,3);

% 
% figure
% set(gcf,'position',[201   201  360  360])
% colormap gray
% subplot('position',[0 0 1 1])
% imagesc(sta(:,:,4));
% set(gca, 'visible', 'off', 'dataaspectratio',[1 1 1])


%%%%%%%%%% PLOT PARASOL CELL STA %%%%%%%%%%
figsize = 135;
figure
set(gcf,'units', 'points','position',[400   400  figsize, figsize])
colormap gray
subplot('position',[0 0 1 1])
imagesc(tmp)
axis([60 176 90 206])
set(gca, 'visible', 'off', 'dataaspectratio',[1 1 1])



%% OFF midget

datarunID = 64;

my_sta=zeros(size(inputs_v,1),sta_params.length);

spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh_v)); % spikes in frames
spikes(spikes<sta_params.length-sta_params.offset)=[];

nspikes = length(spikes);

while ~isempty(spikes)
    [~, ia, ~] = unique(spikes);
    for j=1:sta_params.length
        my_sta(:,sta_params.length-j+1) = my_sta(:, sta_params.length-j+1) +...
            sum(inputs_v(:,spikes(ia)-sta_params.length+j+sta_params.offset),2);
    end
    spikes(ia)=[];
end
my_sta=my_sta/nspikes;

% get voronoi map sta
tmp_map = vormap;
tt=0:max(tmp_map(:));
vorsta=zeros(320,320,sta_params.length);
cnt = 1;
for i=1:1517
    [a, b] = find(tmp_map==tt(i+1));
    if ~isempty(a)
        for j = 1:length(a)
            vorsta(a(j),b(j),:) = my_sta(cnt,:);
        end
    end
    cnt=cnt+1;
end

% figure
% for i=1:6    
%     subplot(2,3,i)
%     colormap gray
%     imagesc(vorsta(:,:,i));
% end
% 
% figure
% vorsta_tmp = vorsta(:, :, 2);
% colormap gray
% imagesc(vorsta_tmp)

tmp = my_sta;
% figure
% plot(my_sta')

thr = robust_std(tmp(:,3))*2;
thr1 = robust_std(tmp(:,2))*2;

tc = mean(tmp(tmp(:,2)<-thr1 & tmp(:,3)>thr,:));

% figure
% plot(tc)
tmp = my_sta(:,1:6).*repmat(tc(:,1:6), size(my_sta,1),1);
tmp = sum(tmp,2);

% get voronoi map sta
tmp_map = vormap;
tt=0:max(tmp_map(:));
vorsta=zeros(320,320);
cnt = 1;
for i=1:1517
    [a, b] = find(tmp_map==tt(i+1));
    if ~isempty(a)
        for j = 1:length(a)
            vorsta(a(j),b(j)) = tmp(cnt);
        end
    end
    cnt=cnt+1;
end

%%%%%%%%%% PLOT MIDGET CELL VORONOI %%%%%%%%%%
figure
set(gcf,'position',[201   201  360  360])
colormap gray
subplot('position',[0 0 1 1])
imagesc(-vorsta)
set(gca, 'visible', 'off', 'dataaspectratio',[1 1 1])
axis([37 176 72 211])


% STA 

sta = squeeze(double(datarun.stas.stas{datarunID}));
% 
% figure
% for i=1:6    
%     subplot(2,3,i)
%     colormap gray
%     imagesc(sta(:,:,i));
% end

tmp = sta(:,:,4);
thr = robust_std(tmp(:))*2;
tmp1 = sta(:,:,3);
thr1 = robust_std(tmp1(:))*2;
coords =  find(tmp<-thr & tmp1>thr1);

tc = [];
for i=1:5
    tmp = sta(:,:,i);
    tc(i)=mean(tmp(coords));
end
% figure
% plot(tc)
tt = shiftdim(repmat(tc',1,320,320),1);
tmp = sta(:,:,1:5).*tt;
tmp = sum(tmp,3);
% 
% figure
% colormap gray
% imagesc(sta(:,:,4));

%%%%%%%%%% PLOT MIDGET CELL STA WITH CONE CENTERS AND VORONOI ON TOP %%%%%%%%%%
figsize = 135;
figure
set(gcf,'units', 'points','position',[400   400  figsize, figsize])
colormap gray
subplot('position',[0 0 1 1])
imagesc(-tmp)
hold on
voronoi_plot_no_centers(datarun1.cones.centers(:,1), datarun1.cones.centers(:,2)) % set colors in the function
set(gca, 'visible', 'off', 'dataaspectratio',[1 1 1])
axis([85 108 133 156])




pp = tmp(94:102,143:148);
pp = pp(pp>(max(pp(:))/2));
tt(1) = mean(pp);

pp = vorsta_tmp(94:102,143:148);
pp = pp(pp>(max(pp(:))/2));
tt(2) = mean(pp);
figure
bar(tt)