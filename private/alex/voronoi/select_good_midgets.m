datarun = load_data('/Volumes/Analysis/2010-09-24-1/d05-36-norefit/data006-from-d05_36/real_sta/data006-from-d05_36');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);


% clean the cells, remove low SNR
sta_all = zeros(320,320);

for cell_type = 1:5
    
    sta_tmp = [];
    cnt = 1; clear rat
    for visionID = datarun.cell_types{cell_type}.cell_ids
        datarunID = find(datarun.cell_ids == visionID);
        sta = squeeze(datarun.stas.stas{datarunID}(:,:,1,:));
        sta = sta(:,:,4);
        if datarun.stas.polarities{datarunID}>0
            sta = -sta;
        end
        [a,b] = find(sta == min(sta(:)),1);
        if ~isempty(a) && a>border && a<320-border && b>border && b<320-border
            cone_template = cone_template + sta(a-border:a+border, b-border:b+border);
        end
        tmin = min(sta(:));
        tmax = max(sta(:));
        rat(cnt) = robust_std(sta(:));%tmin/tmax;
        sta_tmp(:,:,cnt) = sta;
        cnt = cnt+1;
    end
    my_s = find(rat<nanmean(rat)+2*nanstd(rat));
    if cell_type ==1
        sta_all(:,:,end:end+length(my_s)-1)=sta_tmp(:,:,my_s);
    else
        sta_all(:,:,end+1:end+length(my_s))=sta_tmp(:,:,my_s);
    end
    
end

sta_tmp = sta_all;

border = 5;
cone_dist = 1; % min cone distance (everything in between will be wiped out). 1 is 3 stixels.
bb= zeros(size(sta_tmp,3),2,20);
% use both ON and OFF midgets or also nice parasols for this calculation!
i = 1; flag = 1;
while flag
    cone_template = 0;

    cnt = 1; clear rat
    tt = [];
    for j = 1:size(sta_tmp,3)
        sta = sta_tmp(:,:,j);
        [a,b] = find(sta == min(sta(:)),1);
        if ~isempty(a) && a>border && a<320-border && b>border && b<320-border
            cone_template = cone_template + sta(a-border:a+border, b-border:b+border);
            sta_tmp(a-cone_dist:a+cone_dist, b-cone_dist:b+cone_dist,j)= 0;
            tt = [tt; a,b];
        else
            tt = [tt; NaN, NaN];
        end
        
    end
    
    % check for cones which are too far away from the rest for this cell!
    % Stop searching for this cell then.
    
%     tmp = pdist(tt);
%     figure
%     plot(sort(tmp))
    
    bb(:,:,i) = tt;
    
    if i>1
        a = pdist2(bb(:,:,i-1), bb(:,:,i));
        a = diag(a);
        b = a;
        b(isnan(a))=[];
        thresh = 20;%robust_mean(b)+3*robust_std(b);
%         figure
%         plot(a,'*')
%         hold on
%         line([0 120], [thresh thresh])
        tt = find(a>thresh);
        if ~isempty(tt)
            bb(tt,:,i) = NaN;
            sta_tmp(:,:,tt) = NaN; % exclude these cells
        end
        tt = find(a<thresh);
        if isempty(tt)
            flag = 0;
        end
    end
    i = i+1;
%     
%     figure
%     imagesc(cone_template)
%     title(int2str(i))
%     plot(tt(:,1),tt(:,2),'*')
    
end


figure
axis ij
axis([1 320 1 320])
hold on
for j = 1:size(sta_tmp,3)
    my_coord = [];
    for i = 1:size(bb,3)
        my_coord = [my_coord; bb(j,[2,1],i)];
    end
    plot(my_coord(:,1),my_coord(:,2),'-x')
    plot(my_coord(end,1),my_coord(end,2),'-*')
    
end

figure
hold on
for j = 1:10%length(my_cells)
    my_coord = [];
    for i = 1:size(bb,3)
        my_coord = [my_coord; bb(j,[2,1],i)];
    end
    
    figure
    colormap gray
    imagesc(sta_all(:,:,j));
    hold on
    plot(my_coord(:,1),my_coord(:,2),'xr')
    axis([min(my_coord(:,1))-10  max(my_coord(:,1))+10 min(my_coord(:,2))-10  max(my_coord(:,2))+10])
end

a = permute(bb, [1,3,2]);
a = reshape(a, 114*34,2);
a(isnan(a(:,1)),:) = [];

tmp = squareform(pdist(a)) + tril(nan(length(a)));
[c,d] = find(tmp==0);
a(d,:) = [];
tmp = squareform(pdist(a)) + tril(nan(length(a)));


tmp1 = pdist2(a,a);
tmp1 = sort(tmp1);
tmp1(1,:) = [];
td=nanmean(tmp1(1:6,:))

figure
plot(td)
figure
plot(tmp1(1,:))


% check for holes in the cone mosaic! Measure mean distance from each
% pixel to 6 nearest cones. Plot distribution. If it is above certain value
% a cone might be missing. Pull in cells with nearby cones, and look
% specifically at this position to find a potential cone (individual STA
% plots or mean of STAs for cells in this region. Beware of surrounds! -
% take an absolute value of all STAs here?)


% to create a voronoi map with rings:
% 1. Identify cone centers (see above).
% 2. Plot a certain shape (1, cross-like 5, or square 3x3 pixels) around
% the center. this is the main voronoi region.
% 3. Plot a ring around the main Voronoi region, with given width (1,2
% pixels... depends on cone spacing)
% 4. Fill in remaining space with single pixel noise (or more. Depends
% on cone spacing).
% 5. On steps 2 and 3, check that the regions do not overlap: choose
% midpoint and make both regions smaller if too close. If no space for
% rings, leave it alone?





figure
hold on
for i=1:length(my_cells)
    tmp = sta_all(:,:,i);
    [tmp, ic] = sort(tmp(:));
    
    
    plot(tmp(1:200))    
end



figure
colormap gray
imagesc(sta)

figure
colormap gray
imagesc(cone_template)

figure
plot(cone_template(11,:))
hold on
plot(cone_template(:,11))

cone_dist = 3;






datarun1 = load_data('/Volumes/Analysis/2010-09-24-1/d05-36-norefit/data006-from-d05_36/fake_sta/data006-from-d05_36');
datarun1 = load_params(datarun1,'verbose',1);
datarun1 = load_sta(datarun1);

sta = double(squeeze(datarun1.stas.stas{235}));

figure
colormap gray
imagesc(sta(:,:,4))

false_stixels = 100;

[threshold] = sig_stixels_threshold(sta, false_stixels)


datarun = load_data('/Volumes/Analysis/2010-09-24-1/data000/real_sta/data000');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);

datarun1 = load_data('/Volumes/Analysis/2010-09-24-1/data000/fake_sta/data000');
datarun1 = load_params(datarun1,'verbose',1);
datarun1 = load_sta(datarun1);

sta = double(squeeze(datarun.stas.stas{9}));

figure
colormap gray
imagesc(sta(:,:,2,27))

[sig_stixels] = significant_stixels(sta, 'select', 'thresh', 'thresh', threshold)
