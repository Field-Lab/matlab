clear
close all
%number{1} = [481 811 1531 2103 2161 3288 3816 3886 4280 4937 6331 6511];
% number{2} = [185 306 413 996 1058 1759 1821 2061 3396 3593 3830 4656 5342 5733 6064 6213 7148];

number{1} = [81 800 966 1174 2912;...
            3063 3187 3376 3683 3768;...
            4088 4861 5031 5272 5418;...
            6199 6323 7127 7563 nan;...
                   
    ];
%         number{1} = number{1}';
sta_index =2; % 1 for bw, 2 for rgb
threshold = 0.2;

x_plots = size(number{1},1); % number of plots down
y_plots = size(number{1},2);
number{1} = number{1}';
% y_plots = ceil(length(number)/x_plots);
ha = tight_subplot(x_plots, y_plots, [.01 .01],[.01 .01],[.01 .01]);
date = '2015-09-23-7/data031-cf/edited/';
concatname ='data031-cf';
concatname_str = strrep(concatname, '_', '\_')

suptitle([ date, ' ', concatname_str])
for types = 1:size(number,2)
    fig(types) = figure;
    for i = 1:x_plots*y_plots
        if~isnan( number{types}(i))
            i
            clear size_test
            load(['/Volumes/Lab/Users/crhoades/Jitter/', date, '/', concatname,'/Cell ', num2str(number{types}(i)), '.mat'])
            sta = temp;
            
            if size(size(sta),2) == 3 %% BW
                sta = permute(sta, [2,1,3]);
                [junk,start_index] = max(sum(reshape(sta.^2,[],size(sta,3)),1));
                temp = sta(:,:,start_index)';
                
            else
                
                sta = permute(sta, [2,1,3,4]);
                %     [sig_stixels] = significant_stixels(temp, 'select', 'thresh', 'thresh', 4.25);
                [junk,start_index] = max(sum(reshape(sta.^2,[],size(sta,4)),1));
                
                % normalize STA color
                %sta = norm_image(sta);
                
                
                
                % plot spatial sensitivity
                %figure;image(sta(:,:,:,start_index));
                % axis image
                
                %plot_sta_(sta)
                temp = sta(:,:,sta_index,start_index)';
            end
            
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
            
            if size(size(sta),2) == 3 %% BW
                [contour_polygons{types}{i},extras,params] = rf_contours(abs(sta(:,:,start_index)), (max(abs(center_block(:))))*threshold);
                threshold_used{types}(i) = (max(abs(center_block(:))))*threshold;
                axes(ha(i));
                sta = norm_image(sta);
                
                imagesc(sta(:,:,start_index))
                colormap(gray)
            else
                [contour_polygons{types}{i},extras,params] = rf_contours(abs(sta(:,:,sta_index,start_index)), (max(abs(center_block(:))))*threshold);
                threshold_used{types}(i) = (max(abs(center_block(:))))*threshold;
                axes(ha(i));
                sta = norm_image(sta);
                
                image(sta(:,:,:,start_index))
            end
            
            for j= 1:size(contour_polygons{types}{i}{1},2)
                size_test(j) = size(contour_polygons{types}{i}{1}(j).x,2);
            end
            [~, right_ind] = max(size_test);
            hold on
            
            plot(contour_polygons{types}{i}{1}(right_ind).x, contour_polygons{types}{i}{1}(right_ind).y)
            
            
            set(gca, 'ydir', 'reverse')
            axis equal
            axis off
            title(['Cell ',num2str(number{1}(i))])
            
            figure(fig(types))
            hold on
            plot(contour_polygons{types}{i}{1}(right_ind).x, contour_polygons{types}{i}{1}(right_ind).y, 'Linewidth', 2)
            axis equal;
            axis([1, size(sta,2), 1, size(sta,1)])
            box on
            title([ date, ' ', concatname_str])
            set(gca, 'ydir', 'reverse')
            
            [x,y] = poly2cw(contour_polygons{types}{i}{1}(right_ind).x,contour_polygons{types}{i}{1}(right_ind).y);
            area{types}(i) = polyarea(x,y)*(5/1000)^2;%0.005um/pix;
        else
            axes(ha(i));
            axis off
        end
        
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
