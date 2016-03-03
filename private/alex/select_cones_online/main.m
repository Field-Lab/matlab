
%% load data
clear
close all
clear

global cones cone_regions new_cones new_regions

date = '2016-02-17-4';
run = 'data001';
path2load = ['/Volumes/Analysis-1/',date, '/d00-05-norefit/',run,'/',run];
field_size = [300, 400];
scale = 2;

datarun = load_data(path2load);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta','all','keep_java_sta',true);
datarun = set_polarities(datarun);

cell_indices = get_cell_indices(datarun, {1,2,3,4,5});

all_sta = zeros(field_size(1),field_size(2),length(cell_indices));

cnt=1;
t = zeros(1,6);
pols = ones(1, length(cell_indices));
max_pix = zeros(1, length(cell_indices));
for i=cell_indices
    for j=1:6
        sta=squeeze(datarun.stas.stas{i}(:,:,:,j));
        t(j) = max(abs(sta(:)));
    end
    [~, frame] = max(t);
    sta=datarun.stas.stas{i}(:,:,:,frame);
    if max(sta(:))<max(abs(sta(:))) % OFF cell
        pols(cnt) = -1;
    end
    max_pix(cnt) = max(abs(sta(:)))/robust_std(sta(:));
    all_sta(:,:,cnt)=double(squeeze(sta));
    cnt = cnt+1;
end

[~, ic] = sort(max_pix, 'descend');

all_sta = all_sta(:,:,ic);
all_sta = imresize(all_sta,scale); % perhaps nearest?
pols = pols(ic);

%% pre-finding: auto
prelim_peaks = zeros(6000,3);
threshold = 4;
cnt = 0;
tic
for i=1:length(pols)
    sta = pols(i)*all_sta(:,:,i);
    [~, tmp] = find_cones_auto(sta, threshold);
    if ~isempty(tmp)
        prelim_peaks(cnt+1:cnt+size(tmp,1),1:2) = round(tmp);
        prelim_peaks(cnt+1:cnt+size(tmp,1),3) = i;
        cnt = cnt+size(tmp,1);
    end
end
toc
prelim_peaks(cnt+1:end,:) = [];

figure
plot(prelim_peaks(:,1), prelim_peaks(:,2), 'x')

tmp = pdist2(prelim_peaks(:,1:2), prelim_peaks(:,1:2));

ind = 1;
checked_list = [];
full_list = 1:size(all_sta,3);
coll_sta = {};
cct = 1;
while ind<size(all_sta,3) && length(unique(checked_list))<size(all_sta,3)
    
    while ~isempty(find(checked_list==ind, 1))
        ind = ind+1;
    end
    ind
    checked_list = [checked_list ind];
    
    t = find(prelim_peaks(:,3)==ind);
    acc_ic = [];
    for i=1:length(t)
        [tt, ic]=sort(tmp(t(i),:));
        acc_ic = [acc_ic ic(1:5)];
    end
    acc_ic(acc_ic>=t(1) & acc_ic<=t(end)) = [];
    p = unique(acc_ic);
    clear tt
    for i=1:length(p)
        tt(i) = nnz(acc_ic==p(i));
    end
    
    k = find(tt>2);
    sta = (pols(ind)*all_sta(:,:,ind));
    for i = 1:length(k)
        ind2 = prelim_peaks(p(k(i)),3);
        sta = sta + (pols(ind2)*all_sta(:,:,ind2));
        checked_list = [checked_list ind2];
    end
    checked_list = unique(checked_list);
    coll_sta{cct} = sta;
    cct = cct+1;
end

coll_sta = cell2mat(coll_sta);
coll_sta = reshape(coll_sta, 600, 800, cct-1);

pols_a = ones(cct-1,1);
cone_regions = [];
cones = [];
radius = 3
select_cells(coll_sta, pols_a, radius, 1);



%% initial finding

% load('/Volumes/Analysis-1/2016-02-17-4/stimuli/maps/map_data001_test_ready_info.mat');

cone_regions = [];
cones = [];
radius = 3
select_cells(all_sta, pols, radius, 1);


% eliminate absolute duplicates
flag = 1;
while flag
    a = triu(pdist2(cones, cones));
    a(a==0) = 100;
    [r,c] = find(a<1);
    if isempty(r)
        flag = 0;
    else
        doubled_cone = [];
        for i=1:length(c)
            if isempty(find(r==c(i), 1))
                doubled_cone = [doubled_cone c(i)];
            end
        end
        cones(unique(doubled_cone),:) = [];
        cone_regions(unique(doubled_cone)) = [];
    end
end

path2save=['/Volumes/Analysis-1/', date, '/stimuli/maps/map_', run, '_test_20160302'];
if ~isdir(path2save)
    mkdir(path2save);
end
save([path2save '_info_delete'],'cell_indices', 'ic', 'cones', 'cone_regions')

%% control: repeat until no double cones

% important parameters:
min_overlap = 0.2; % if cones' centers are closer than this distance, one cone is made out of them automatically
max_overlap = 6; % if cone centers are farther than this distance, they will be dubbed separate cones


flag = 1;
p=0;
while flag
    a = triu(pdist2(cones, cones));
    tmp = zeros(length(cones))+100;
    tmp = tril(tmp);
    a = a+tmp;
    [r,c] = find(a<min_overlap);
    if isempty(r)
        flag = 0;
    else
        doubled_cone = [];
        for i=1:length(c)
            if isempty(find(r==c(i), 1))
                doubled_cone = [doubled_cone c(i)];
            end
        end
        p = p+length(doubled_cone);
        cones(unique(doubled_cone),:) = [];
        cone_regions(unique(doubled_cone)) = [];
    end
end
p

% [a, b] = find(cones(:,1)>400)
% cones(1781,:) = [];
% cone_regions(1781) = [];

figure
hold on
axis ij
for i=1:length(cone_regions)
    plot(cone_regions{i}(:,1), cone_regions{i}(:,2));
end


w = zeros(size(cone_regions,2),size(all_sta,3));
for i=1:size(cone_regions,2)
    t = round(cones(i,:));
    cr = cone_regions{i};
    [X, Y] = meshgrid(linspace(t(1)-10,t(1)+10, 21), linspace(t(2)-10,t(2)+10, 21));
    [in, on] = inpolygon(X, Y, cr(:,1), cr(:,2));
    in = find(in);
    tmp = 0;
    for j=1:length(in)
        tmp = tmp + squeeze(all_sta(Y(in(j)), X(in(j)), :));
    end
    w(i,:) = tmp/length(in);
end

for i=1:size(all_sta,3)
    w(:,i) = w(:,i)*pols(i);
end

tmp_w = w./repmat(max(w), size(w,1),1);


a = triu(pdist2(cones, cones));
a(a==0) = 100;
ct = zeros(size(a,1),2);
ct(:,1) = 1:size(a,1);
for i=1:size(a,1)
    t = [find(a(i,:)<max_overlap) i];
    if any(ct(t,2))
        p = [];
        for j=t
            if  ct(j,2)~=0
                p = [p; ct(ct(:,2) == ct(j,2))];
            else
                p = [p; j];
            end
        end
        p = unique(p);
        val = sort(ct(p,2));
        val = val(find(val,1));
        ct(p,2) = val;
    elseif length(t)>1
        ct(t,2) = max(ct(:,2))+1;
    end
    p=unique(ct(:,2));
    clear c
    for k=1:length(p)
        c(k) = length(find(ct(:,2)==p(k)));
    end
    if any(c==1)
        c;
        disp(i)
    end
end

cone_lists = {};
cnt = 1;
for i=1:max(ct(:,2))
    p = find(ct(:,2)==i);
    if ~isempty(p)
        cone_lists{cnt} = p;
        cnt = cnt+1;
    end
end

double_cones = cell2mat(cone_lists');

length(double_cones)

radius = 3;
weight_threshold = 0.25;
correct_cones(tmp_w, all_sta, pols, cones, cone_regions, cone_lists, radius, weight_threshold)

cones(double_cones,:) = [];
cone_regions(double_cones) = [];
cones = [cones; new_cones];
cone_regions = [cone_regions new_regions];


%% construct map
map = zeros(600,800);
for i=1:size(cone_regions,2)
    t = round(cones(i,:));
    cr = cone_regions{i};
    [X, Y] = meshgrid(linspace(t(1)-10,t(1)+10, 21), linspace(t(2)-10,t(2)+10, 21));
    [in, on] = inpolygon(X, Y, cr(:,1), cr(:,2));
    in = find(in);
    for j=1:length(in)
        map(Y(in(j)), X(in(j))) = i;
    end
end
figure
imagesc(map)
tmp = map;
tmp(tmp>1) = 1;
figure
imagesc(tmp)


%% save map
path2save=['/Volumes/Analysis-1/', date, '/stimuli/maps/map_', run, '_test_20160302'];
dlmwrite([path2save '.txt'], map, 'delimiter', '\t', 'newline', 'pc');
save([path2save '_info'],'cell_indices', 'ic', 'cones', 'cone_regions')

savedMap=dlmread([path2save '.txt']);
figure
imagesc(savedMap)


%% control!

for i=1:10
    tmp = pols(i)*all_sta(:,:,i);
    
    comb = zeros(600,800,3);
    comb(:,:,1) = map;
    comb(comb>0) = 0.5;
    comb(:,:,2) = tmp/max(tmp(:));
    figure
    imagesc(comb)
end



%% track cones



marks = zeros(size(all_sta,1),size(all_sta,2),size(all_sta,4));
threshold = 3.8;
for ind=1:size(all_sta,4)
    
    if rgb
        sta=all_sta(:,:,:,ind);
        sta2show = 1/max(abs(sta(:)));
        sta2show = sta*sta2show/2+0.5;
    else
        sta=all_sta(:,:,ind);
        sta2show = sta;
    end
    [cons, cone_peaks] = find_cones_auto(pols(ind)*sta, threshold);
    for i=1:size(cons,1)
        marks(cons(i,1),cons(i,2),ind) = 1;
    end
end


radius = 1.5;

% do automatic adjustment??? before manual clicking??
track_cones(all_sta, marks, pols, cones, cone_regions, radius*2, 17)

tic
[new_cones, new_regions] = track_cones_auto(all_sta, marks, pols, cones, radius*2);
toc


figure
hold on
axis ij
for i=1:length(new_regions)
    if ~isempty(new_regions{i})
        plot(new_regions{i}(:,1), new_regions{i}(:,2));
    end
end
figure
hold on
axis ij
for i=1:length(cone_regions)
    plot(cone_regions{i}(:,1), cone_regions{i}(:,2), 'k');
end

figure
plot(new_cones(:,1)/2, new_cones(:,2)/2, 'xr');
hold on
plot(cones(:,1), cones(:,2), '+k');



% construct map
radius = 3;
x = 0;
y = 0;
th = 0:pi/250:2*pi;
xunit = round(((radius * cos(th) + x)*2))/2;
yunit = round(((radius * sin(th) + y)*2))/2;
a = [xunit; yunit]';
m = [];
for j=1:499
    if xunit(j)==xunit(j+1) && yunit(j)==yunit(j+1)
        m = [m j];
    end
end
a(m,:) = [];

[X, Y] = meshgrid(linspace(x-10,x+10, 21), linspace(y-10,y+10, 21));
in = inpolygon(X, Y, a(:,1), a(:,2));
[rows, cols] = find(in);
rows = rows-11;
cols = cols-11;
figure
plot(a(:,1), a(:,2))

map = zeros(600,800);
for i=1:size(new_regions,2)
    if ~isempty(new_regions{i})
        
        x = new_cones(i,1);
        y = new_cones(i,2);
        
        trows = rows+x;
        tcols = cols+y;
        linearInd = sub2ind([600,800], tcols, trows);
        map(linearInd) =i;
    end
    
end

figure
imagesc(map)


for i=1:20
    tmp = pols(i)*all_sta(:,:,2,i);
    tmp =  imresize(tmp,scale,'method', 'nearest')/max(tmp(:));
    [j, k] = find(tmp == max(tmp(:)),1);
    
    comb = zeros(600,800,3);
    comb(:,:,1) = map;
    comb(comb>0) = 0.5;
    comb(:,:,2) = tmp;
    figure
    imagesc(comb)
    axis([k-20 k+20 j-20 j+20])
end











check_missed =[];
for i=1:length(prelim_peaks)
    t = prelim_peaks(i,1:2);
    if map(t(2), t(1))==0
        check_missed = [check_missed prelim_peaks(i,3)];
    end
end
check_missed = unique(check_missed);



for i=check_missed(11:30)
    tmp = pols(i)*all_sta(:,:,i);
    
    comb = zeros(600,800,3);
    comb(:,:,1) = map;
    comb(comb>0) = 0.5;
    comb(:,:,2) = tmp/max(tmp(:));
    figure
    imagesc(comb)
end

