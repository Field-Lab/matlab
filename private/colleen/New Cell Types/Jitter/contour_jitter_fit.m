%% Plot STA overlaid with largest contour found at a user defined threshold

clear
close all



% vision cell ids in an array which is how they show up in subplots 
% Use nan's if you have a non rectangular number of cells ie [1 2; 3 nan]
% if testing multiple cell types, use number{1} = , number{2} =  etc. 
number{1} = [527 902 1568;...
            2371 3186 4041]
        
        
date = '2015-09-23-7/data031-cf/edited/';
concatname ='data031-cf';

sta_index =2; % 1 for bw, 2 for rgb

% Set manually
threshold = 0.2;

%% End of inputs
x_plots = size(number{1},1); % number of plots down
y_plots = size(number{1},2);
number{1} = number{1}';

ha = tight_subplot(x_plots, y_plots, [.01 .01],[.01 .01],[.01 .01]);

concatname_str = strrep(concatname, '_', '\_');

suptitle([ date, ' ', concatname_str])
for types = 1:size(number,2)
    fig(types) = figure;
    for i = 1:x_plots*y_plots
        if~isnan( number{types}(i))
            disp(i)
            clear size_test
            
            % path to the sta
            % loads as the variable 'temp'
            load(['/Volumes/Lab/Users/crhoades/Jitter/', date, '/', concatname,'/Cell ', num2str(number{types}(i)), '.mat'])
            sta = temp;
            
            if size(size(sta),2) == 3 %% BW
                sta = permute(sta, [2,1,3]);
                [junk,start_index] = max(sum(reshape(sta.^2,[],size(sta,3)),1));
                temp = sta(:,:,start_index)';
                
            else %% RGB
                
                sta = permute(sta, [2,1,3,4]);
                [junk,start_index] = max(sum(reshape(sta.^2,[],size(sta,4)),1));
                temp = sta(:,:,sta_index,start_index)';
            end
            
            % Find a box around the center
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
%             noise_pixels = temp([1:left, right:end], [1:bottom, top:end]);
%             noise_level{types}(i) = median(noise_pixels(:));
            
            if size(size(sta),2) == 3 %% BW
                % normalize threshold by the max amplitude in the center
                [contour_polygons{types}{i},extras,params] = rf_contours(abs(sta(:,:,start_index)), (max(abs(center_block(:))))*threshold);
                threshold_used{types}(i) = (max(abs(center_block(:))))*threshold;
                axes(ha(i));
                sta = norm_image(sta);
                
                imagesc(sta(:,:,start_index))
                colormap(gray)
            else %% RGB
                 % normalize threshold by the max amplitude in the center

                [contour_polygons{types}{i},extras,params] = rf_contours(abs(sta(:,:,sta_index,start_index)), (max(abs(center_block(:))))*threshold);
                threshold_used{types}(i) = (max(abs(center_block(:))))*threshold;
                axes(ha(i));
                sta = norm_image(sta);
                
                image(sta(:,:,:,start_index))
            end
            
            % Find the largest contour
            for j= 1:size(contour_polygons{types}{i}{1},2)
                size_test(j) = size(contour_polygons{types}{i}{1}(j).x,2);
            end
            [~, right_ind] = max(size_test);
            hold on
            
            % plot the largest contour
            plot(contour_polygons{types}{i}{1}(right_ind).x, contour_polygons{types}{i}{1}(right_ind).y)
            
            
            set(gca, 'ydir', 'reverse')
            axis equal
            axis off
            title(['Cell ',num2str(number{types}(i))])
            
            % Summary figure of all the contours in the population of cells
            figure(fig(types))
            hold on
            plot(contour_polygons{types}{i}{1}(right_ind).x, contour_polygons{types}{i}{1}(right_ind).y, 'Linewidth', 2)
            axis equal;
            axis([1, size(sta,2), 1, size(sta,1)])
            box on
            title([ date, ' ', concatname_str])
            set(gca, 'ydir', 'reverse')
            
            [x,y] = poly2cw(contour_polygons{types}{i}{1}(right_ind).x,contour_polygons{types}{i}{1}(right_ind).y);
            % area of the contour
            final_area{types}(i) = polyarea(contour_polygons{types}{i}{1}(right_ind).x, contour_polygons{types}{i}{1}(right_ind).y);

        else %nan in vision id array
            axes(ha(i));
            axis off
        end
        
    end
end


%% If you want to plot RF area of two populations of cells 
% figure;
% hold on
% [nb,xb] =hist(area{1});
% createPatches(xb,nb,[1 1 1],0.005,0.05);
% [nb,xb] =hist(area{2});
% createPatches(xb,nb,[0 0 0],0.005, 0.75);


