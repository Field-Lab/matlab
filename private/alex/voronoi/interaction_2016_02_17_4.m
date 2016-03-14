%% load stuff

local_path = '/Volumes/Analysis-1/';

vormap = load([local_path, '2016-02-17-4/stimuli/maps/map_data001_from_data000_600_1500.txt']);

figure
colormap gray
imagesc(vormap)

temp = uint8(zeros(size(vormap)));
voronoi_contours = cell(max(vormap(:)),2);
figure
colormap gray
td = vormap;
td(vormap>0) = 1;
imagesc(td)
hold on
for i=1:max(vormap(:))
    tmpmap=temp;
    if ~isempty(find(vormap==i,1))
        tmpmap(vormap==i)=1;
        dd = imresize(tmpmap,5,'method', 'nearest');
        [r, c] = find(dd,1);
        contour = bwtraceboundary(dd,[r c],'W',8,Inf,'counterclockwise');
        contour= round(contour/5)+0.5;
        plot(contour(:,2), contour(:,1), 'linewidth', 2)
        voronoi_contours{i, 1} = contour(:,2);
        voronoi_contours{i, 2} = contour(:,1);
    end
end

%Coarse run run 1
datarunc = load_data([local_path, '2016-02-17-4/d00-05-norefit/data000/data000']);
datarunc = load_params(datarunc,'verbose',1);
datarunc = load_sta(datarunc);


%SC run 1
datarun = load_data([local_path, '2016-02-17-4/d00-05-norefit/data001/data001']);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);


%SC run OTHER DATE
datarun7 = load_data([local_path, '2016-02-17-7/d01-13-norefit/data002/data002']);
datarun7 = load_params(datarun7,'verbose',1);
datarun7 = load_sta(datarun7);


%SC run 1 FALSE
datarunf = load_data([local_path, '2016-02-17-4/d00-05-norefit/data001_false/data001']);
datarunf = load_params(datarunf,'verbose',1);
datarunf = load_sta(datarunf);

%SC run 2
datarun2 = load_data([local_path, '2016-02-17-4/d00-05-norefit/data005/data005']);
datarun2 = load_params(datarun2,'verbose',1);
datarun2 = load_sta(datarun2);

% voronoi run
vorrun = load_data([local_path, '2016-02-17-4/d00-05-norefit/data004/data004']);
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_sta(vorrun);
vorrun = load_neurons(vorrun);
[inputs, refresh, duration] = get_wn_movie_ath(vorrun, 'BW-1-1-0.48-11111-1340x1-60.35.xml');

%% cone profiles: before
all_sta = zeros(600,800,length(datarun.cell_ids));
cnt = 1;
pols = zeros(1,length(datarun.cell_ids));
for i = 1:length(datarun.cell_ids)
    tmp = datarun.stas.stas{i};
    t = zeros(1,6);
    for j=1:6
        sta=squeeze(tmp(:,:,:,j));
        t(j) = max(abs(sta(:)));
    end
    [~, frame] = max(t);
    sta=tmp(:,:,:,frame);
    if max(sta(:))<max(abs(sta(:))) % OFF cell
        pols(cnt) = -1;
    end
    all_sta(:,:,cnt)=imresize(double(squeeze(sta)), 2, 'method', 'nearest');
    cnt = cnt+1;
end

cone_portrait = cell(max(vormap(:)),2);
for i=1:max(vormap(:))
    
    [a,b] = find(vormap==i);
    tmp = 0;
    for j=1:length(a)
        tmp = tmp+squeeze(all_sta(a(j),b(j),:));
    end
    tmp = tmp'.*pols;
    cone_portrait{i,1} = find(tmp>(max(tmp)*0.7));
    cone_portrait{i,2} = [a b];
end

figure
bord = 15;
cnt = 1;
for i=28:36
    tmp = cone_portrait{i,1};
    if ~isempty(tmp)
        tt = 0;
        for j=1:length(tmp)
            tt = tt+pols(tmp(j))*all_sta(:,:,tmp(j));
        end
        subplot(3,3,cnt)
        imagesc(tt)
        tmp =  cone_portrait{i,2};
        axis([min(tmp(:,2))-bord  max(tmp(:,2))+bord min(tmp(:,1))-bord  max(tmp(:,1))+bord ])
        hold on
        plot(voronoi_contours{i, 1}, voronoi_contours{i, 2}, 'color', 'r')
        title(int2str(j))
    end
    cnt = cnt+1;
end


bord = 2;
cnt = 1;
combo = zeros(600,800);
for i=1:max(vormap(:))
    tmp = cone_portrait{i,1};
    if ~isempty(tmp)
        tt = 0;
        for j=1:length(tmp)
            tt = tt+pols(tmp(j))*all_sta(:,:,tmp(j));
        end
        tmp =  cone_portrait{i,2};
        combo(min(tmp(:,1))-bord:max(tmp(:,1))+bord, min(tmp(:,2))-bord:max(tmp(:,2))+bord) = tt(min(tmp(:,1))-bord:max(tmp(:,1))+bord, min(tmp(:,2))-bord:max(tmp(:,2))+bord);
    end
    cnt = cnt+1;
end
figure
imagesc(combo)
hold on
for i=1:max(vormap(:))
    plot(voronoi_contours{i, 1}, voronoi_contours{i, 2}, 'color', 'r')
end

%% cone profiles: after
good_cells = [541, 556, 571, 647, 796, 1276, 1367, 1486, 1816, 1892, 1921, 2087, 2117,...
    2161, 2192, 2222, 2386, 2658, 3001, 3031, 3046, 3182, 3258, 3272, 3721, 4111, 4486,...
    4966, 4967, 5058, 5341, 5613, 5822, 5914, ...
    6017, 6061, 6244, 6256, 6332, 6483, 6572, 6991, 7156, 7246, 7502, 7652, 7667, 7697];
all_sta = zeros(600,800,length(good_cells));
cnt = 1;
pols = zeros(1,length(good_cells));
max_pics = pols;
for i = 1:length(datarun.cell_ids)
    if ~isempty(find(good_cells == datarun.cell_ids(i),1))
        tmp = datarun.stas.stas{i};
        t = zeros(1,6);
        for j=1:6
            sta=squeeze(tmp(:,:,:,j));
            t(j) = max(abs(sta(:)));
        end
        [~, frame] = max(t);
        sta=tmp(:,:,:,frame);
        if max(sta(:))<max(abs(sta(:))) % OFF cell
            pols(cnt) = -1;
        end
        all_sta(:,:,cnt)=imresize(double(squeeze(sta)), 2, 'method', 'nearest');
        tmp = sta*pols(cnt);
        max_pics(cnt) = max(tmp(:));
        cnt = cnt+1;
    end
end

cone_portrait = cell(max(vormap(:)),2);
for i=1:max(vormap(:))
    
    [a,b] = find(vormap==i);
    tmp = 0;
    for j=1:length(a)
        tmp = tmp+squeeze(all_sta(a(j),b(j),:));
    end
    tmp = tmp'.*pols;
    cone_portrait{i,1} = find(tmp>(max(tmp)*0.85) & tmp>max_pics*0.9);
    cone_portrait{i,2} = [a b];
end


bord = 2;
cnt = 1;
combo = zeros(600,800);
for i=1:max(vormap(:))
    tmp = cone_portrait{i,1};
    if ~isempty(tmp)
        tt = 0;
        for j=1:length(tmp)
            tt = tt+pols(tmp(j))*all_sta(:,:,tmp(j));
        end
        tmp =  cone_portrait{i,2};
        combo(min(tmp(:,1))-bord:max(tmp(:,1))+bord, min(tmp(:,2))-bord:max(tmp(:,2))+bord) = tt(min(tmp(:,1))-bord:max(tmp(:,1))+bord, min(tmp(:,2))-bord:max(tmp(:,2))+bord);
    end
    cnt = cnt+1;
end
figure
imagesc(combo)
hold on
for i=1:max(vormap(:))
    plot(voronoi_contours{i, 1}, voronoi_contours{i, 2}, 'color', 'r')
end


for k = 1:9:100
    
    figure
    set(gcf, 'position', [-1580         125        1085         966])
    bord = 15;
    cnt = 1;
    for i=k:k+8
        tmp = cone_portrait{i,1};
        if ~isempty(tmp)
            tt = 0;
            for j=1:length(tmp)
                tt = tt+pols(tmp(j))*all_sta(:,:,tmp(j));
            end
            subplot(3,3,cnt)
            imagesc(tt)
            tmp =  cone_portrait{i,2};
            axis([min(tmp(:,2))-bord  max(tmp(:,2))+bord min(tmp(:,1))-bord  max(tmp(:,1))+bord ])
            hold on
            plot(voronoi_contours{i, 1}, voronoi_contours{i, 2}, 'color', 'r')
            title(int2str(j))
        end
        cnt = cnt+1;
    end
    
end



%% cone profiles: after
all_sta = zeros(600,800,length(datarun.cell_ids));
cnt = 1;
pols = zeros(1,length(datarun.cell_ids));
max_pics = pols;
for i = 1:length(datarun.cell_ids)
    tmp = datarun.stas.stas{i};
    t = zeros(1,6);
    for j=1:6
        sta=squeeze(tmp(:,:,:,j));
        t(j) = max(abs(sta(:)));
    end
    [~, frame] = max(t);
    sta=tmp(:,:,:,frame);
    if max(sta(:))<max(abs(sta(:))) % OFF cell
        pols(cnt) = -1;
    end
    all_sta(:,:,cnt)=imresize(double(squeeze(sta)), 2, 'method', 'nearest');
    tmp = sta*pols(cnt);
    max_pics(cnt) = max(tmp(:));
    cnt = cnt+1;
end

cone_portrait = cell(max(vormap(:)),2);
for i=1:max(vormap(:))
    
    [a,b] = find(vormap==i);
    tmp = 0;
    for j=1:length(a)
        tmp = tmp+squeeze(all_sta(a(j),b(j),:));
    end
    tmp = tmp'.*pols;
    cone_portrait{i,1} = find(tmp>(max(tmp)*0.85) & tmp>max_pics*0.9);
    cone_portrait{i,2} = [a b];
end



bord = 2;
cnt = 1;
combo = zeros(600,800);
for i=1:max(vormap(:))
    tmp = cone_portrait{i,1};
    if ~isempty(tmp)
        tt = 0;
        for j=1:length(tmp)
            tt = tt+pols(tmp(j))*all_sta(:,:,tmp(j));
        end
        tmp =  cone_portrait{i,2};
        combo(min(tmp(:,1))-bord:max(tmp(:,1))+bord, min(tmp(:,2))-bord:max(tmp(:,2))+bord) = tt(min(tmp(:,1))-bord:max(tmp(:,1))+bord, min(tmp(:,2))-bord:max(tmp(:,2))+bord);
    end
    cnt = cnt+1;
end
figure
imagesc(combo)
hold on
for i=1:max(vormap(:))
    plot(voronoi_contours{i, 1}, voronoi_contours{i, 2}, 'color', 'r')
end

%% cone profiles: after
all_sta = zeros(300,400,length(datarun.cell_ids));
cnt = 1;
pols = zeros(1,length(datarun.cell_ids));
max_pics = pols;
for i = 1:length(datarun.cell_ids)
    tmp = datarun.stas.stas{i};
    t = zeros(1,6);
    for j=1:6
        sta=squeeze(tmp(:,:,:,j));
        t(j) = max(abs(sta(:)));
    end
    [~, frame] = max(t);
    sta=tmp(:,:,:,frame);
    if max(sta(:))<max(abs(sta(:))) % OFF cell
        pols(cnt) = -1;
    end
    all_sta(:,:,cnt)=double(squeeze(sta));
%     all_sta(:,:,cnt)=imresize(double(squeeze(sta)), 2, 'method', 'nearest');
    tmp = sta*pols(cnt);
    max_pics(cnt) = max(tmp(:));
    cnt = cnt+1;
end

all_sta1 = zeros(300,400,length(datarun7.cell_ids));
cnt = 1;
pols = zeros(1,length(datarun7.cell_ids));
max_pics = pols;
for i = 1:length(datarun7.cell_ids)
    tmp = datarun7.stas.stas{i};
    t = zeros(1,6);
    for j=1:6
        sta=squeeze(tmp(:,:,:,j));
        t(j) = max(abs(sta(:)));
    end
    [~, frame] = max(t);
    sta=tmp(:,:,:,frame);
    if max(sta(:))<max(abs(sta(:))) % OFF cell
        pols(cnt) = -1;
    end
    all_sta1(:,:,cnt)=double(squeeze(sta));
%     all_sta(:,:,cnt)=imresize(double(squeeze(sta)), 2, 'method', 'nearest');
    tmp = sta*pols(cnt);
    max_pics(cnt) = max(tmp(:));
    cnt = cnt+1;
end


%% sta from coarse
all_staC = zeros(24,32,length(datarun.cell_ids));
cnt = 1;
pols = zeros(1,length(datarun.cell_ids));
max_pics = pols;
for i = 1:length(datarun.cell_ids)
    tmp = datarunc.stas.stas{i};
    t = zeros(1,30);
    for j=1:30
        sta=squeeze(tmp(:,:,2,j));
        t(j) = max(abs(sta(:)));
    end
    [~, frame] = max(t);
    sta=tmp(:,:,2,frame);
    if max(sta(:))<max(abs(sta(:))) % OFF cell
        pols(cnt) = -1;
    end
    all_staC(:,:,cnt)=double(squeeze(sta));
%     all_sta(:,:,cnt)=imresize(double(squeeze(sta)), 2, 'method', 'nearest');
    tmp = sta*pols(cnt);
    max_pics(cnt) = max(tmp(:));
    cnt = cnt+1;
end




tic
comb_sta = zeros(600,800);
for i=1:600
    for k=1:800
        tmp = squeeze(all_sta(i,k,:))'.*pols;
        tmp = find(tmp>robust_std(tmp)*10 & tmp>max_pics*0.2);
        
        if ~isempty(tmp)
            comb_sta(i,k) = mean(pols(tmp).*squeeze(all_sta(i,k,tmp))');
        end
        
    end
end
toc
figure
imagesc(comb_sta)

a = get_cell_indices(datarun, {4});

tmp = sum(all_sta(:,:,a),3);
figure
imagesc(tmp)

tmp = sum(all_sta(:,:,a(1)),3);
figure
imagesc(tmp)

tmp = all_sta(1:100,1:100,a(1));
mean(tmp(:))
b = mean(all_sta(:,:,a(1)));
plot(b)
figure
b = mean(all_sta(:,:,a(1))');
plot(b)

figure
hold on
b = 0;
for i=1:100
    b = squeeze(all_sta(1,i,a))+b;
end
plot(sort(b)/100)


b = sort(squeeze(all_sta1(64,313,a)));
sum(b)
hold on
plot(b, 'k', 'linewidth', 3)
b = sort(squeeze(all_sta(90,304,a)));
sum(b)
plot(b, 'b', 'linewidth', 3)
b = sort(squeeze(all_sta1(73,318,a)));
sum(b)
plot(b, 'g', 'linewidth', 3)
b = sort(squeeze(all_sta1(1,2,a)));
sum(b)
plot(b, 'r', 'linewidth', 3)
b = sort(squeeze(all_sta(end,end,a)));
sum(b)
plot(b, 'm', 'linewidth', 3)

t = [];
for i=1:55
    t(i) = sum(b(i:end-i+1));
end
% figure
plot(t, 'g', 'linewidth', 3)
hold on

tmp = sum(all_sta(:,:,a),3);

figure
hist(tmp(:),100)

tmp = sum(all_sta(1:20,1:20,a),3);
figure
hist(tmp(:),20)

figure
cnt = 1;
for i=1:5
    for j=1:5
        tmp = squeeze(all_sta(i,j,a));
        [p ic ]= sort(tmp);
        subplot(5,5,cnt)
        plot(-abs(tmp(ic)))
        line([55 55], [min(-abs(tmp(ic))) 0], 'color', 'r')
        cnt = cnt+1;
    end
end


i = [313 313 303 304 318 310 320 326 324 317]
k = [64 63 90 90 73 94 97 69 74 87]


figure
cnt = 1;
for j = 1:length(i)-1
    tmp = squeeze(all_sta(k(j),i(j),a));
    [p ic ]= sort(tmp);
    subplot(3,3,cnt)
    plot(-abs(tmp(ic(3:end))))
    line([54 54], [min(-abs(tmp(ic(3:end)))) 0], 'color', 'r')
    cnt = cnt+1;
end


t = [];
cnt = 1;
for i=[1:10, 291:300]
    for j=[1:10, 350:400]
        tmp = squeeze(all_sta(i,j,a));        
        [p ic ]= sort(tmp);
        t(cnt)=find(p>0,1);
        cnt = cnt+1;
    end
end
figure
hist(t, 50)

figure
p = [];
for i=1:100
    c = randperm(110);
    b = sort(squeeze(all_sta1(1,2,a(c(1:100)))));
    p(i) = mean(b);
end
plot(p)
hold on

a = get_cell_indices(datarun, {4});

figure
subplot(2,2,2)
p = [];
for i=1:400
    t = double(rand(45000*length(a),1)>0.5);
    t(t==0) = -1;
    t(t>0) = 1;
    p(i) = mean(t);
end
plot(p)
axis([1 400 -1e-2 1e-2]) 
title('simulation, [-1 1]')

subplot(2,2,2)
p = [];
for i=1:400
    t = double(rand(45000*length(a),1)>0.5);
    t(t==0) = -0.48;
    t(t>0) = 0.48;
    p(i) = mean(t);
end
plot(p)
axis([1 400 -1e-2 1e-2]) 
title('simulation, [-0.48 0.48]')

subplot(2,2,3)
p = [];cnt = 1;
for i=[1:10, 291:300]
    for j=[1:10, 391:400]
        b = squeeze(all_sta1(i,j,a));
        p(cnt) = mean(b);
        cnt = cnt+1;
    end
end
plot(p)
axis([1 400 -1e-2 1e-2]) 
title('data, true seed')

figure
imagesc(all_sta(:,:,1))

subplot(2,2,4)
k = [];cnt = 1;
for i=[1:10, 291:300]
    for j=[1:10, 391:400]
        b = squeeze(all_sta(i,j,a));
        k(cnt) = mean(b);
        cnt = cnt+1;
    end
end
plot(p)
axis([1 400 -1e-2 1e-2]) 
title('data, false seed')

corr(k',p')
a = get_cell_indices(datarun,{4});
tmp = all_sta(:,:,a(1));
figure;imagesc(tmp)
td = tmp([1:10, 291:300], [1:10, 391:400]);

 b = all_sta([1:10, 291:300],[1:10, 391:400],get_cell_indices(datarun,{4}));
 b = reshape(b, 400,110);
plot(mean(b'))

b2 = mean(b1,3)
b2 = (mean(b1,3)+mean(b,3))/2;
 tmp([1:10, 291:300], [1:10, 391:400]) =  tmp([1:10, 291:300], [1:10, 391:400])-b2;

 b1 = all_sta1([1:10, 291:300],[1:10, 391:400],get_cell_indices(datarun7,{4}));
 b1 = reshape(b1, 400,27);
plot(mean(b1'))

a= tmp([1:10, 291:300], [1:10, 391:400]);
mean(a(:))
std(a(:))


% cell 1
cell_inds = get_cell_indices(datarun,{4});
tmp = all_sta(:,:,cell_inds(1));
tmp2 = tmp;
b = all_sta(50:110,300:340,cell_inds(2:end));
b = mean(b,3);
tmp2(50:110,300:340) =  tmp(50:110,300:340)-b;

figure;
imagesc(tmp)

% cell 2
cell_inds = get_cell_indices(datarun,{4});
tmp = all_sta(:,:,cell_inds(2));
tmp2 = tmp;
b = all_sta(:,:,cell_inds([1, 3:end]));
b = mean(b,3);
tmp2 =  tmp-b;

figure;
imagesc(tmp)

% combine 2 cells
cell_inds = get_cell_indices(datarun,{4});
tmp = all_sta(:,:,cell_inds(1));
tmp1 = all_sta(:,:,cell_inds(2));
b = all_sta(:,:,cell_inds(2:end));
c = all_sta(:,:,cell_inds([1, 3:end]));
b = mean(b,3);
c = mean(c,3);
tmp2 =  tmp+tmp1-b-c;
tmp3 = tmp + tmp1;
figure;
imagesc(tmp2)
title('subtract')
figure;
imagesc(tmp3)
title('raw')

% combine all
tmp1 = 0;
tmp2 = 0;
for i=1:10
    cell_inds = get_cell_indices(datarun,{4});
    tmp = all_sta(:,:,cell_inds(i));
    cell_inds(i) = [];
    b = all_sta(:,:,cell_inds);
    b = mean(b,3);
    tmp1 = tmp1 + tmp - b;
    tmp2 = tmp2 + tmp;
end
figure;
imagesc(tmp1)
title('subtract')
figure;
imagesc(tmp2)
title('raw')

% combine all
tmp1 = 0;
tmp2 = 0;
cell_inds = get_cell_indices(datarun,{4});
tmp2 = sum(all_sta(:,:,cell_inds(1:10)),3);

for i=1:10
    tmp = all_sta(:,:,cell_inds(i));
    cell_inds(i) = [];
    b = all_sta(:,:,cell_inds);
    b = mean(b,3);    
end
figure;
imagesc(tmp1)
title('subtract')
figure;
imagesc(tmp2)
title('raw')

% combine all - non-cone regions
% noise map
cell_inds = get_cell_indices(datarun7,{2,3,4});
noise_sta = 0;
for i=1:77 
    tmp = all_sta1(:,:,cell_inds(i));
    [a,b] = find(abs(tmp)==max(abs(tmp(:))));
%     if a>41 & a<258 & b>41 & b<358
% %         tmp(a-40:a+40,b-40:b+40) = 0;
        noise_sta = noise_sta+tmp;    
%     else
%         i
%     end
    
end
noise_sta = noise_sta/77;
figure;
imagesc(noise_sta)

cell_inds = get_cell_indices(datarun7,{4});
noise_sta_offm = 0;
for i=1:27 
    tmp = all_sta1(:,:,cell_inds(i));
    [a,b] = find(abs(tmp)==max(abs(tmp(:))));
%     tmp(a-40:a+40,b-40:b+40) = 0;
    noise_sta_offm = noise_sta_offm+tmp;    
end
noise_sta_offm = noise_sta_offm/27;
figure;
imagesc(noise_sta_offm)


cell_inds = get_cell_indices(datarun7,{4});
tmp = sum(all_sta1(:,:,cell_inds(1:5)),3);
figure
imagesc(tmp)
tmp = tmp- noise_sta_offm*5;
figure
imagesc(tmp)

title('subtract')
figure;
imagesc(tmp)
title('raw')


cell_inds = get_cell_indices(datarun7,{2});
for i=1:10
    tmp = sum(all_sta1(:,:,cell_inds(i)),3);
    tmp1 = tmp-noise_sta;
    tmp2 = tmp - noise_sta_offm;
    figure
    set(gcf, 'position', [-1870         550        1871         521])
    ax1 = subplot(1,3,1);
    colormap gray
    imagesc(tmp1)
    title(['denoised with all , ', num2str(robust_mean(tmp1(:))), ' + ', num2str(robust_std(tmp1(:)))])
    ax2 = subplot(1,3,2);
    colormap gray
    imagesc(tmp2)
    title(['denoised with off midgets, ', num2str(robust_mean(tmp1(:))), ' + ', num2str(robust_std(tmp2(:)))])
    ax3 = subplot(1,3,3);
    colormap gray
    imagesc(tmp)
    title(['raw, ', num2str(robust_mean(tmp1(:))), ' + ', num2str(robust_std(tmp(:)))])
    linkaxes([ax1,ax2,ax3],'xy')
end


% combine all - COARSE
% noise map
cell_inds = get_cell_indices(datarun,{1,2,3,4});
noise_sta = 0;
for i=1:208 
    tmp = all_staC(:,:,cell_inds(i));
    [a,b] = find(abs(tmp)==max(abs(tmp(:))));
%     if a>41 & a<258 & b>41 & b<358
% %         tmp(a-40:a+40,b-40:b+40) = 0;
        noise_sta = noise_sta+tmp;    
%     else
%         i
%     end
    
end
noise_sta = noise_sta/208;
figure;
imagesc(noise_sta)

cell_inds = get_cell_indices(datarun,{1,2,3,4});
noise_sta_emp = 0;
cnt = 208;
for i=1:208 
    tmp = all_staC(:,:,cell_inds(i));
    [a,b] = find(abs(tmp)==max(abs(tmp(:))));
    if a>4 & a<20 & b>4 & b<28
        tmp(a-4:a+4,b-4:b+4) = 0;
        noise_sta_emp = noise_sta_emp+tmp;
    else
        cnt = cnt-1;
    end
    
end
noise_sta_emp = noise_sta_emp/191;
figure;
imagesc(noise_sta_emp)



cell_inds = get_cell_indices(datarunc,{1});
for i=1:10
    tmp = sum(all_staC(:,:,cell_inds(i)),3);
    tmp1 = tmp-noise_sta;
    tmp2 = tmp - noise_sta_emp;
    figure
    set(gcf, 'position', [-1870         550        1871         521])
    ax1 = subplot(1,3,1);
    colormap gray
    imagesc(tmp1)
    title(['denoised with all , ', num2str(robust_mean(tmp1(:))), ' + ', num2str(robust_std(tmp1(:)))])
    ax2 = subplot(1,3,2);
    colormap gray
    imagesc(tmp2)
    title(['denoised with all emp, ', num2str(robust_mean(tmp1(:))), ' + ', num2str(robust_std(tmp2(:)))])
    ax3 = subplot(1,3,3);
    colormap gray
    imagesc(tmp)
    title(['raw, ', num2str(robust_mean(tmp1(:))), ' + ', num2str(robust_std(tmp(:)))])
    linkaxes([ax1,ax2,ax3],'xy')
end



figure
subplot(2,2,4)
p = [];cnt = 1;
for i=151:170
    for j=200:219
        b = squeeze(all_sta(i,j,a));
        p(cnt) = mean(b);
        cnt = cnt+1;
    end
end
plot(p)
axis([1 400 -1e-2 1e-2]) 
title('data, false seed')



state = Init_RNG_JavaStyle(11111);

p = [];
for j=1:400
    rand_seq = [];
    for i=1:21532*110
        rand_seq(i) = random_uint16(state)>32768;
    end
    t = double(rand_seq);
    t(t==0) = -0.48;
    t(t>0) = 0.48;
    p(j) = mean(t);
end
figure
plot(p)
axis([1 400 -1e-2 1e-2])
title('sim, [-0.48 0.48], java rand')



figure
hold on
for k = 1:4
    
    a = get_cell_indices(datarun, {k});
    p = [];cnt = 1;
    for i=[1:10, 291:300]
        for j=[1:10, 391:400]
            b = squeeze(all_sta(i,j,a));
            p(cnt) = mean(b);
            cnt = cnt+1;
        end
    end
    plot(p)
end

figure
hold on
for k = 1:4
    
    a = get_cell_indices(datarun, {k});
    p = [];cnt = 1;
    for i=1:300
        for j=1:400
            b = squeeze(all_sta(i,j,a));
            p(cnt) = mean(b);
            cnt = cnt+1;
        end
    end
    plot(p)
end

figure
hold on
for k = 1:4
    
    a = get_cell_indices(datarun, {k});
    p = [];cnt = 1;
    for i=1:24
        for j=1:32
            b = squeeze(all_staC(i,j,a));
            p(cnt) = mean(b);
            cnt = cnt+1;
        end
    end
    plot(p)
end




figure
hold on
for k = 1:4
    
    a = get_cell_indices(datarun, {k});
    p = [];cnt = 1;
    for i=[1:10, 291:300]
        for j=[1:10, 391:400]
            i
            j
            b = squeeze(all_sta(i,j,a(1:5)));
            p(cnt) = mean(b);
            cnt = cnt+1;
        end
    end
    plot(p)
end
legend('ONp', 'OFFp', 'ONm', 'OFFm')


figure
hold on

a = get_cell_indices(datarun, {4});
p = [];cnt = 1;
for i=[1:10, 291:300]
    for j=[1:10, 391:400]
        b = squeeze(all_sta(i,j,a(1:5)));
        p(cnt) = mean(b);
        cnt = cnt+1;
    end
end
plot(p)
legend('ONp', 'OFFp', 'ONm', 'OFFm')




a = get_cell_indices(datarun, {4});
c = get_cell_indices(datarun, {3});
p = [];cnt = 1;
for i=[1:10, 291:300]
    for j=[1:10, 391:400]
        b = squeeze(all_sta(i,j,[a c]));
        p(:,cnt) = b;
        cnt = cnt+1;
    end
end
p = p';
tmp = corr(p);
tmp = triu(tmp,1);
t = tmp(1:110,111:end);
k = tmp(1:110,1:110);
k(k==0) = [];
figure
imagesc(t)
mean(t(:))
std(t(:))

mean(k(:))
std(k(:))

figure
plot(sort(t))

figure
imagesc(all_sta(:,:,c(end)))

    

figure
hold on
for k=1:2:10
    for j=1:2:10
        b = sort(squeeze(all_sta(k,j,a)));
        t = [];
        for i=1:55
            t(i) = sum(b(i:end-i+1-10));
        end
        plot(t)
    end
end



figure
hold on
for k=60:2:70
    for j=309:2:319
        b = sort(squeeze(all_sta(k,j,a)));
        t = [];
        for i=1:55
            t(i) = sum(b(i:end-i+1-10));
        end
        plot(t)
    end
end




tmp = all_sta;
figure
hold on
for i=1:100
    b = sort(squeeze(tmp(i,1,a)));
    plot(b)
end
[b, ic1] = sort(squeeze(tmp(126,626,a)));
plot(b, 'k', 'linewidth', 3)
[b, ic2] = sort(squeeze(tmp(173,634,a)));
plot(b, 'r', 'linewidth', 3)
b = squeeze(tmp(126,626,a));
plot(-b(ic1), 'g', 'linewidth', 3)

bord = 2;
cnt = 1;
combo = zeros(600,800);
for i=1:max(vormap(:))
    tmp = cone_portrait{i,1};
    if ~isempty(tmp)
        tt = 0;
        for j=1:length(tmp)
            tt = tt+pols(tmp(j))*all_sta(:,:,tmp(j));
        end
        tmp =  cone_portrait{i,2};
        combo(min(tmp(:,1))-bord:max(tmp(:,1))+bord, min(tmp(:,2))-bord:max(tmp(:,2))+bord) = tt(min(tmp(:,1))-bord:max(tmp(:,1))+bord, min(tmp(:,2))-bord:max(tmp(:,2))+bord);
    end
    cnt = cnt+1;
end
figure
imagesc(combo)
hold on
for i=1:max(vormap(:))
    plot(voronoi_contours{i, 1}, voronoi_contours{i, 2}, 'color', 'r')
end






%% plot stuff
save_path = [local_path, '2016-02-17-4/data004_d00-05/'];
nbins_cone1 = 8;
nbins_cone2 = 8;
sta_params.length = 20;
sta_params.offset = 0;
fraction = 0.9;


for kkk= vorrun.cell_types{1, 1}.cell_ids(1:end)  %vorrun.cell_ids(22:end)
    close all
    datarunID = find(vorrun.cell_ids==kkk);
    
    for i=1:length(vorrun.cell_types)
        if ~isempty(find(vorrun.cell_types{i}.cell_ids==kkk, 1))
            cell_type = vorrun.cell_types{i}.name;
            break
        end
    end
    
    visionID = vorrun.cell_ids(datarunID);
    raw_sta = squeeze(vorrun.stas.stas{datarunID});
    %     figure;
    %     plot(raw_sta')
    thresh = mean([robust_std(raw_sta(:,5)),robust_std(raw_sta(:,10))])*5;
    
    a = min(raw_sta(:));
    b = max(raw_sta(:));
    if abs(a)>abs(b) % OFF cell
        pol = -1;
        [~,pos] = find(raw_sta==a,1);
    else % ON cell
        pol = 1;
        [~,pos] = find(raw_sta==b,1);
        a = b;
    end
    tmp = raw_sta * pol;
    
    spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh)); % spikes in frames
    spikes(spikes<sta_params.length-sta_params.offset)=[];
    
    if abs(a) > thresh*1.5 && length(spikes)>300
        
        [~, cones] = expand_voronoi_sta(raw_sta, vormap);
        
        center_cones = find(sum(tmp(:,pos-1:pos+1)>thresh,2));
        center_cones(isnan(cones(center_cones,1)))=[];
        
        if length(center_cones)>1
            
            
            %     figure
            %     plot(raw_sta(center_cones,:)')
            timecourse = mean(tmp(center_cones,10:end));
            %     figure
            %     plot(timecourse)
            
            conv_sta = sum(tmp(:,10:end).*repmat(timecourse,size(tmp,1),1),2);
            if length(center_cones)>30
                [~, a] = sort(conv_sta);
                center_cones = sort(a(end-29:end));
            end
            timecourse = mean(tmp(center_cones,10:end));
            conv_sta = sum(tmp(:,10:end).*repmat(timecourse,size(tmp,1),1),2);
            
            % prepare voronoi and single cone STAs
            [full_sta, cones] = expand_voronoi_sta(conv_sta, vormap);
            full_sta = pol*full_sta/(max(full_sta(:))*1.5);
            sta1 = squeeze(datarun.stas.stas{datarunID});
            sta1 = imresize(sta1(:,:,4), 2, 'method', 'nearest');
            sta2 = squeeze(datarun2.stas.stas{datarunID});
            sta2 = imresize(sta2(:,:,4), 2, 'method', 'nearest');
            
            
            
           
            [unbiased_sta, gensig_bins, nonlinearity]=unbiased_STA(inputs(center_cones,1:50000), spikes, fraction, sta_params);
            
            %         figure
            %         plot(unbiased_sta')
            
            filt_inputs = zeros(length(center_cones), size(inputs,2)-sta_params.length+1);
            cnt = 1;
            for current_cone=center_cones'
                filt_inputs(cnt,:)=conv(inputs(current_cone,:), unbiased_sta(cnt,:),'valid');
                cnt=cnt+1;
            end
            spikes_tmp = spikes;
            spikes_tmp(spikes<sta_params.length) = [];
            
            offline_contrast_response(filt_inputs, spikes_tmp-sta_params.length+1,...
                nbins_cone1, nbins_cone2, center_cones, vormap, cones, full_sta, datarunID,...
                kkk, [save_path int2str(nbins_cone2) '_bins/', cell_type '/'], sta1, sta2, voronoi_contours, pol);
        end
    end
end



