clear
close all
tic
fixed_threshold = 0.11;
starting_threshold = 0.08;
threshold_increment = 0.005;
area_min = 100; % contours smaller than this area (in pixels squared) won't be considered
size_of_null = 300; % number of examples to make up the null distribution
mode = 'unfixed_thresh'; % fixed_thres... or something else
alpha = 1e-60; % for t test

% number must be an array; Fill in nan where cells aren't needed
number{1} = [527; 902; 1622;...
    ];

number{2} = [4594; 5193 ;5493;...
    ];


sta_index =2; % 1 for bw, 2 for rgb (green channel)

x_plots{1} = size(number{1},1); % number of plots down
y_plots{1} = size(number{1},2);
x_plots{2} = size(number{2},1); % number of plots down
y_plots{2} = size(number{2},2);
number{1} = number{1}';
number{2} = number{2}';


date = '2016-04-21-1/data006-cf/edited';
concatname ='data006-cf';
concatname_str = strrep(concatname, '_', '\_')

suptitle([ date, ' ', concatname_str])
for types = 1:size(number,2)
    fig(types) = figure;
    fig(types+2) = figure;
    figure;
    ha = tight_subplot(x_plots{types}, y_plots{types}, [.01 .01],[.01 .01],[.01 .01]);
    suptitle([ date, ' ', concatname_str])
    
    for i = 1:x_plots{types}*y_plots{types} % interate over the number of cells
        if~isnan( number{types}(i)) % if that part of the plot has a cell in it
            i
            clear size_test
            
            % load the jitter STA
            load(['/Volumes/Lab/Users/crhoades/Jitter_backup/', date, '/', concatname,'/Cell ', num2str(number{types}(i)), '.mat'])
            sta = temp;
            % Normalize by the max (in the future, divided by the number of
            % spikes if that wasn't done when the STA was computed
            sta = sta./max(abs(sta(:)));
            
            h = 0; % hypothesis true or false
            p = 0; % p value vector
            
            
            if size(size(sta),2) == 3 %% BW
                sta = permute(sta, [2,1,3]); % switch x and y dimensions
                [junk,start_index] = max(sum(reshape(sta.^2,[],size(sta,3)),1)); % find the peak frame
                temp = sta(:,:,start_index)';
                
            else
                
                sta = permute(sta, [2,1,3,4]); % switch x and y dimensions
                [junk,start_index] = max(sum(reshape(sta.^2,[],size(sta,4)),1)); % find the peak frame
                temp = sta(:,:,sta_index,start_index)';
            end
            
            
            norm_sta = norm_image(sta);
            
            if strcmp(mode, 'fixed_thresh')
                
                % put in new plot handles
                
                
                threshold = fixed_threshold;
                
                if size(size(sta),2) == 3 %% BW
                    % Find the contours at the given threshold level
                    [contour_polygons{types+2}{i},extras,params] = rf_contours(abs(sta(:,:,start_index)), threshold);
                    
                    cla(ha(i));
                    axes(ha(i));
                    
                    hold off
                    
                    
                    
                    imagesc(norm_sta(:,:,start_index))
                    colormap(gray)
                else % RGB
                    % Find the contours at the given threshold level
                    [contour_polygons{types+2}{i},extras,params] = rf_contours(abs(sta(:,:,sta_index,start_index)), threshold);
                    
                    cla(ha(i));
                    axes(ha(i));
                    
                    hold off
                    
                    
                    imagesc(norm_sta(:,:,:,start_index))
                end
                hold on
                axis equal
                axis off
                final_area{types}(i) = 0;
                
                for j = 1:size(contour_polygons{types+2}{i}{1},2)
                    area(j) = polyarea(contour_polygons{types+2}{i}{1}(j).x, contour_polygons{types+2}{i}{1}(j).y);
                    if area(j) < area_min % eliminate contours that are too small
                        contour_polygons{types+2}{i}{1}(j).x = nan;
                        contour_polygons{types+2}{i}{1}(j).y = nan;
                    end
                    if ~isnan(contour_polygons{types+2}{i}{1}(j).x)
                        
                        inside = 0;
                        % test if there are contours within other contours
                        for k = 1:size(contour_polygons{types+2}{i}{1},2) % increment through all the real contours
                            if ~isnan(contour_polygons{types+2}{i}{1}(j).x)
                                if k ~=j % if the contour to test isn't the same as the contour of comparison
                                    % inside will be greater than zero if the
                                    % contour of interest is a small contour
                                    % inside a bigger real one
                                    inside = inside + sum(inpolygon(contour_polygons{types+2}{i}{1}(j).x, contour_polygons{types+2}{i}{1}(j).y, contour_polygons{types+2}{i}{1}(k).x, contour_polygons{types+2}{i}{1}(k).y));
                                end
                            end
                            
                        end
                        if inside == 0
                            % only plot the outermost contour
                            
                            
                            axes(ha(i))
                            
                            
                            
                            plot(contour_polygons{types+2}{i}{1}(j).x, contour_polygons{types+2}{i}{1}(j).y, 'linewidth',2)
                            title(['Cell ',num2str(number{types}(i))])
                            
                            figure(fig(types+2))
                            hold on
                            plot(contour_polygons{types+2}{i}{1}(j).x, contour_polygons{types+2}{i}{1}(j).y, 'linewidth',2)
                            final_area{types}(i) = final_area{types}(i) + polyarea(contour_polygons{types+2}{i}{1}(j).x, contour_polygons{types+2}{i}{1}(j).y);
                            
                        end
                        
                    end
                end
                
                
                
                axis equal;
                axis([1, size(sta,2), 1, size(sta,1)])
                box on
                title([ date, ' ', concatname_str])
                set(gca, 'ydir', 'reverse')
                
                
                
                
                
                
                
                
            else % Trying to find the optimal threshold for all the cells
                
                threshold = starting_threshold;
                
                while length(h(~isnan(h))) ~= sum(h(~isnan(h))) % while the number of "true" hypotheses isn't all of them
                    clear h
                    clear p
                    % raise the threshold
                    threshold = threshold + threshold_increment;
                    
                    if size(size(sta),2) == 3 %% BW
                        % Find the contours at the given threshold level
                        [contour_polygons{types}{i},extras,params] = rf_contours(abs(sta(:,:,start_index)), threshold);
                        
                        
                        cla(ha(i));
                        axes(ha(i))
                        hold off
                        
                        imagesc(norm_sta(:,:,start_index))
                        colormap(gray)
                    else
                        % Find the contours at the given threshold level
                        
                        [contour_polygons{types}{i},extras,params] = rf_contours(abs(sta(:,:,sta_index,start_index)), threshold);
                        cla(ha(i));
                        axes(ha(i))
                        
                        image(norm_sta(:,:,:,start_index))
                    end
                    
                    BW = zeros(size(sta,1), size(sta,2),size(contour_polygons{types}{i}{1},2)); % used to make a mask to compute power inside the contour
                    
                    % this loop computes the power inside all the contours
                    % found
                    for j= 1:size(contour_polygons{types}{i}{1},2) % interate through number of contours
                        % compute the area of each contour
                        area(j) = polyarea(contour_polygons{types}{i}{1}(j).x, contour_polygons{types}{i}{1}(j).y);
                        if area(j) < area_min % eliminate contours that are too small
                            contour_polygons{types}{i}{1}(j).x = nan;
                            contour_polygons{types}{i}{1}(j).y = nan;
                            size_test(j) = 0;
                            BW(:,:,j) = 0;
                            power(j) = 0;
                        else
                            
                            size_test(j) = size(contour_polygons{types}{i}{1}(j).x,2);
                            
                            % Make a mask that is just what's inside the
                            % contour
                            BW(:,:,j) = poly2mask(contour_polygons{types}{i}{1}(j).x,contour_polygons{types}{i}{1}(j).y, size(sta,1), size(sta,2));
                            % multiply the mask (1s and 0s) with the sta
                            % and divide by the number of pixels inside the
                            % contour to find the power inside the contour
                            power(j) = sum(sum(BW(:,:,j).*sta(:,:,sta_index, start_index)))./sum(sum(BW(:,:,j)));
                        end
                        
                    end
                    
                    % This look decides if a contour is a real cell
                    for j = 1:size(contour_polygons{types}{i}{1},2)
                        % if the contour was too small, ignore it
                        if isnan(power(j)) | power(j) == 0
                            h(j) = nan;
                            p(j) = nan;
                        else
                            
                            % Draw a bunch of random locations of the
                            % contour
                            rand_seed = randi(size(sta,1)*size(sta,2), size_of_null*10,1); % draw more random numbers than I need in case the null contour overlaps with the real one
                            counter = 1;
                            z = 0;
                            % Keep going until I run out of random
                            % locations or I have size_of_null number of
                            % null example powers
                            while counter < size_of_null+1 && z+1<length(rand_seed)
                                z = z+1; % counter for random draws
                                [y,x] = ind2sub([size(sta,1), size(sta,2)],rand_seed(z)); %y < 320 x<640
                                diff_x = contour_polygons{types}{i}{1}(j).x(1)-x; %.x goes to 640
                                diff_y = contour_polygons{types}{i}{1}(j).y(1)-y; %.y goes to 320
                                
                                % Find a random location for the contour
                                placement_x = contour_polygons{types}{i}{1}(j).x - diff_x;
                                placement_y = contour_polygons{types}{i}{1}(j).y - diff_y;
                                
                                if sum(placement_x < 1 | placement_x > size(sta,2)) == 0 && sum(placement_y < 1 | placement_y > size(sta,1))==0 && sum(inpolygon(placement_x, placement_y, contour_polygons{types}{i}{1}(j).x, contour_polygons{types}{i}{1}(j).y)) == 0 % make sure the contour doesn't exceed the size of the sta or overlap the real location of the contour
                                    moved_mask = poly2mask(placement_x,placement_y, size(sta,1), size(sta,2)); % put the mask over the fake contour location
                                    null_power(j,counter) = sum(sum(moved_mask.*sta(:,:,sta_index, start_index)))./sum(moved_mask(:)); %compute the power in the fake contour
                                    counter = counter+1; % increment the number of successful members of null distribution
                                    
                                    
                                    
                                end
                            end
                            % perform a t test to see if the real contour
                            % is sufficiently outside the null distribution
                            [h(j),p(j)] = ttest2(null_power(j,:), power(j), 'Alpha', alpha);
                        end
                        
                        
                    end
                    
                    % Display the p values of those contours that were
                    % determined to be real cells
                    right_ind = find(p <=alpha);
                    disp(p(right_ind))
                    
                end
                
                % now that we reached the threshold where all the contours found were
                % "real", plot them
                hold on
                for j = 1:length(right_ind)
                    
                    if length(right_ind) == 1 %only one real contour
                        plot(contour_polygons{types}{i}{1}(right_ind(j)).x, contour_polygons{types}{i}{1}(right_ind(j)).y, 'linewidth',2)
                    else
                        inside = 0;
                        % test if there are contours within other contours
                        for k = 1:length(right_ind) % increment through all the real contours
                            if k ~=j % if the contour to test isn't the same as the contour of comparison
                                % inside will be greater than zero if the
                                % contour of interest is a small contour
                                % inside a bigger real one
                                inside = inside + sum(inpolygon(contour_polygons{types}{i}{1}(right_ind(j)).x, contour_polygons{types}{i}{1}(right_ind(j)).y, contour_polygons{types}{i}{1}(right_ind(k)).x, contour_polygons{types}{i}{1}(right_ind(k)).y));
                            end
                        end
                        if inside == 0
                            % only plot the outermost contour
                            plot(contour_polygons{types}{i}{1}(right_ind(j)).x, contour_polygons{types}{i}{1}(right_ind(j)).y, 'linewidth',2)
                            
                        end
                    end
                    
                    hold on
                    axis equal
                    axis off
                end
                
                
                % store the threshold optimized for each cell
                threshold_used{types}(i) = threshold;
                
                set(gca, 'ydir', 'reverse')
                axis equal
                axis off
                title(['Cell ',num2str(number{1}(i))])
                
                % Produce second figure with all the contours of a cell
                % type overlayed
                figure(fig(types))
                hold on
                
                for j = 1:length(right_ind)
                    if length(right_ind) == 1
                        % only one real contour
                        plot(contour_polygons{types}{i}{1}(right_ind(j)).x, contour_polygons{types}{i}{1}(right_ind(j)).y, 'linewidth',2)
                        % compute the area for each cell
                        final_area{types}(i) = polyarea(contour_polygons{types}{i}{1}(right_ind(j)).x, contour_polygons{types}{i}{1}(right_ind(j)).y);
                        
                    else
                        % test if there are contours within other contours
                        
                        inside = 0;
                        for k = 1:length(right_ind)
                            if k ~=j % if the contour to test isn't the same as the contour of comparison
                                % inside will be greater than zero if the
                                % contour of interest is a small contour
                                % inside a bigger real one
                                
                                inside = inside + sum(inpolygon(contour_polygons{types}{i}{1}(right_ind(j)).x, contour_polygons{types}{i}{1}(right_ind(j)).y, contour_polygons{types}{i}{1}(right_ind(k)).x, contour_polygons{types}{i}{1}(right_ind(k)).y));
                            end
                        end
                        final_area{types}(i) = 0;
                        if inside == 0
                            % only plot the outermost contour
                            
                            plot(contour_polygons{types}{i}{1}(right_ind(j)).x, contour_polygons{types}{i}{1}(right_ind(j)).y, 'linewidth',2)
                            final_area{types}(i) = final_area{types}(i) + polyarea(contour_polygons{types}{i}{1}(right_ind(j)).x, contour_polygons{types}{i}{1}(right_ind(j)).y);
                            
                        end
                    end
                    
                    hold on
                    axis equal
                    axis off
                end
                
                
                axis equal;
                axis([1, size(sta,2), 1, size(sta,1)])
                box on
                title([ date, ' ', concatname_str])
                set(gca, 'ydir', 'reverse')
                
            end
            
        else % There is no cell id given for a particular location
            axes(ha(i));
            axis off
        end
        
    end
end

% If not fixed threshold mode, compute the area of all the cells again,
% using the max threshold.
min_x = nan;
min_y = nan;
max_x = nan;
max_y = nan;
if ~strcmp(mode, 'fixed_thresh')
    threshold = max([threshold_used{1}, threshold_used{2}]);
    for types = 1:size(number,2)
        if exist('ha')
            hb= ha;
        end
        
        figure;
        ha = tight_subplot(x_plots{types}, y_plots{types}, [.01 .01],[.01 .01],[.01 .01]);
        suptitle([ date, ' ', concatname_str, ' threshold: ', num2str(threshold)])
        for i = 1:x_plots{types}*y_plots{types} % interate over the number of cells
            if~isnan(number{types}(i)) % if that part of the plot has a cell in it
                i
                clear size_test
                
                % load the jitter STA
                load(['/Volumes/Lab/Users/crhoades/Jitter_backup/', date, '/', concatname,'/Cell ', num2str(number{types}(i)), '.mat'])
                sta = temp;
                % Normalize by the max (in the future, divided by the number of
                % spikes if that wasn't done when the STA was computed
                sta = sta./max(abs(sta(:)));
                
                h = 0; % hypothesis true or false
                p = 0; % p value vector
                
                
                if size(size(sta),2) == 3 %% BW
                    sta = permute(sta, [2,1,3]); % switch x and y dimensions
                    [junk,start_index] = max(sum(reshape(sta.^2,[],size(sta,3)),1)); % find the peak frame
                    temp = sta(:,:,start_index)';
                    
                else
                    
                    sta = permute(sta, [2,1,3,4]); % switch x and y dimensions
                    [junk,start_index] = max(sum(reshape(sta.^2,[],size(sta,4)),1)); % find the peak frame
                    temp = sta(:,:,sta_index,start_index)';
                end
                
                
                norm_sta = norm_image(sta);
                
                
                % put in new plot handles
                
                
                
                if size(size(sta),2) == 3 %% BW
                    % Find the contours at the given threshold level
                    [contour_polygons{types+2}{i},extras,params] = rf_contours(abs(sta(:,:,start_index)), threshold);
                    
                    cla(ha(i));
                    axes(ha(i));
                    
                    hold off
                    
                    
                    
                    imagesc(norm_sta(:,:,start_index))
                    colormap(gray)
                else % RGB
                    % Find the contours at the given threshold level
                    [contour_polygons{types+2}{i},extras,params] = rf_contours(abs(sta(:,:,sta_index,start_index)), threshold);
                    
                    cla(ha(i));
                    axes(ha(i));
                    
                    hold off
                    
                    
                    imagesc(norm_sta(:,:,:,start_index))
                end
                hold on
                axis equal
                axis off
                final_area{types}(i) = 0;
                
                for j = 1:size(contour_polygons{types+2}{i}{1},2)
                    area(j) = polyarea(contour_polygons{types+2}{i}{1}(j).x, contour_polygons{types+2}{i}{1}(j).y);
                    if area(j) < area_min % eliminate contours that are too small
                        contour_polygons{types+2}{i}{1}(j).x = nan;
                        contour_polygons{types+2}{i}{1}(j).y = nan;
                    end
                    if ~isnan(contour_polygons{types+2}{i}{1}(j).x)
                        
                        inside = 0;
                        % test if there are contours within other contours
                        for k = 1:size(contour_polygons{types+2}{i}{1},2) % increment through all the real contours
                            if ~isnan(contour_polygons{types+2}{i}{1}(j).x)
                                if k ~=j % if the contour to test isn't the same as the contour of comparison
                                    % inside will be greater than zero if the
                                    % contour of interest is a small contour
                                    % inside a bigger real one
                                    inside = inside + sum(inpolygon(contour_polygons{types+2}{i}{1}(j).x, contour_polygons{types+2}{i}{1}(j).y, contour_polygons{types+2}{i}{1}(k).x, contour_polygons{types+2}{i}{1}(k).y));
                                end
                            end
                            
                        end
                        if inside == 0
                            % only plot the outermost contour
                            
                            
                            axes(ha(i))
                            
                            
                            
                            plot(contour_polygons{types+2}{i}{1}(j).x, contour_polygons{types+2}{i}{1}(j).y, 'linewidth',2)
                            title(['Cell ',num2str(number{types}(i))])
                            
                            figure(fig(types+2))
                            hold on
                            plot(contour_polygons{types+2}{i}{1}(j).x, contour_polygons{types+2}{i}{1}(j).y, 'linewidth',2)
                            final_area{types}(i) = final_area{types}(i) + polyarea(contour_polygons{types+2}{i}{1}(j).x, contour_polygons{types+2}{i}{1}(j).y);
                            min_x = min(min_x, min(contour_polygons{types+2}{i}{1}(j).x));
                            min_y = min(min_y, min(contour_polygons{types+2}{i}{1}(j).y));
                            max_x = max(max_x, max(contour_polygons{types+2}{i}{1}(j).x));
                            max_y = max(max_y, max(contour_polygons{types+2}{i}{1}(j).y));
                            
                        end
                        
                    end
                end
                
                
                axis equal;
                axis([1, size(sta,2), 1, size(sta,1)])
                box on
                title([ date, ' ', concatname_str, ' threshold: ', num2str(threshold)])
                
                set(gca, 'ydir', 'reverse')
                
            else % There is no cell id given for a particular location
                axes(ha(i));
                axis off
            end
            
        end
        
        
    end
    for types = 1:size(number,2)
        
        for i = 1:x_plots{types}*y_plots{types}
            
            if types == 1
                axes(hb(i))
            else
                axes(ha(i))
            end
            if~isnan( number{types}(i)) % if that part of the plot has a cell in it
                
                hold on
                rectangle('Position', [min_x, min_y, max_x - min_x, max_y-min_y])
                
            end
            
            
        end
        
        
        figure(fig(types+2))
        hold on
        rectangle('Position', [min_x, min_y, max_x - min_x, max_y-min_y])
        
        
        
    end
    
end








toc

figure;

plot(final_area{1}*5.5^2/1000^2, ones(size(final_area{1})), 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'w')
hold on;
plot(final_area{2}*5.5^2/1000^2, 2*ones(size(final_area{2})) , 'ok', 'MarkerSize', 5, 'MarkerFaceColor', 'k')
set(gca, 'ytick', [1 2])
set(gca, 'yticklabel', {['ON Parasol: n = ', num2str(length(final_area{1}))]; ['OFF Parasol: n = ', num2str(length(final_area{2}))]});
set(gca, 'ylim', [0 3]);
xlabel('area (mm^2)')
title({[ date, ' ', concatname_str], 'ON and OFF Parasol Cell Area Comparison', ['Threshold: ', num2str(threshold)]})

% set(bh, 'facecolor', [0 0 0 0.25])
%             hold on
%
%             [nb,xb] =hist(final_area{1}*5.5^2/1000^2);
%             createPatches(xb,nb,[1 1 1],0.005,0.05);
%             [nb,xb] =hist(final_area{2}*5.5^2/1000^2);
%             % [bh]= bar(xb,nb);
%             createPatches(xb,nb,[0 0 0],0.005, 0.75);

% datarun = get_rf_contours(datarun, cell_specification, contour_levels, varargin)

