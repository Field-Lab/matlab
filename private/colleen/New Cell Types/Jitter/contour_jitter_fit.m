clear
close all
%number{1} = [481 811 1531 2103 2161 3288 3816 3886 4280 4937 6331 6511];
number{2} = [1936 2011 3166 3901 4113 4981 5012 6526 7051 7456];
number{2} = [ 4113 4981 5012 ];

threshold = 0.2;
for types = 1:size(number,2)
fig(types) = figure;
for i = 1:size(number{types},2)
    i
    clear size_test
    load(['/Volumes/Lab/Users/crhoades/Jitter/2016-02-17-6/data026cf_spikes/Cell ', num2str(number{types}(i)), '.mat'])
    sta = temp;
    
    
    sta = permute(sta, [2,1,3,4]);
    %     [sig_stixels] = significant_stixels(temp, 'select', 'thresh', 'thresh', 4.25);
    [junk,start_index] = max(sum(reshape(sta.^2,[],size(sta,4)),1));
    
    % normalize STA color
    %sta = norm_image(sta);
    
    
    
    % plot spatial sensitivity
    %figure;image(sta(:,:,:,start_index));
    % axis image
    
    %plot_sta_(sta)
    temp = sta(:,:,2,start_index)';
    %   sta = sta./max(temp(:));
    [val,ind] = max(abs(temp(:)));
    [x,y] = ind2sub(size(temp), ind);
    left = x-50;
    right = x+49;
    bottom = y-50;
    top = y+49;
    if x < 50
        left = 0;
    end
    if x > size(temp,1)-49
        right = size(temp,1)-49;
    end
    if y < 50
        bottom = 50;
    end
    if y > size(temp,2)-49
        top =size(temp,2)-49;
    end
    
    center_block = temp(left: right, bottom:top);
    noise_pixels = temp([1:left, right:end], [1:bottom, top:end]);
    noise_level{types}(i) = median(noise_pixels(:));
    
    [contour_polygons{types}{i},extras,params] = rf_contours(abs(sta(:,:,2,start_index)), (max(abs(center_block(:))))*threshold);
    threshold_used{types}(i) = (max(abs(center_block(:))))*threshold;
    sta = norm_image(sta);
    figure; image(sta(:,:,:,start_index))
    for j= 1:size(contour_polygons{types}{i}{1},2)
        size_test(j) = size(contour_polygons{types}{i}{1}(j).x,2);
    end
    [~, right_ind] = max(size_test);
    hold on
    
    plot(contour_polygons{types}{i}{1}(right_ind).x, contour_polygons{types}{i}{1}(right_ind).y)
    
    
    set(gca, 'ydir', 'reverse')
    
    
    figure(fig(types))
    hold on
    plot(contour_polygons{types}{i}{1}(right_ind).x, contour_polygons{types}{i}{1}(right_ind).y)
    axis equal;
    [x,y] = poly2cw(contour_polygons{types}{i}{1}(right_ind).x,contour_polygons{types}{i}{1}(right_ind).y);
    area{types}(i) = polyarea(x,y)*(5/1000)^2;%0.005um/pix;
end
end


figure;

% set(bh, 'facecolor', [0 0 0 0.25])
hold on

[nb,xb] =hist(area{1});
createPatches(xb,nb,[1 1 1],0.005,0.05);
[nb,xb] =hist(area{2});
% [bh]= bar(xb,nb);
createPatches(xb,nb,[0 0 0],0.005, 0.75);

% datarun = get_rf_contours(datarun, cell_specification, contour_levels, varargin)
    
