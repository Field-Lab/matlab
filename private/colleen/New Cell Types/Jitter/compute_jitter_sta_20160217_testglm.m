%% --------------- Compute STA of focal white noise in MATLAB (not Java) -------------------
%% Function: Generate STAs for cells that were stimulated with the Voronoi stimulus

%% How to use:
%     1) Called from sta_focal_WN
%     2) Need to have spikes loaded and computes an STA by parsing the movie file correctly based on the number of cells targeted in the run


%% Potential problems:
%     1) Could have a problem parsing movie file correctly for multiple
%     cells (limited testing)
%     2) Binning problems with matlab STA code



%% Inputs
% datarun : generated in sta_focal_WN by load_data
% mdf_file : white noise xml
% num_frames : integer of max frame number (eg 30)
% spikes: n x1 vector of spike times in seconds. Obtained from datarun
% plotting: 1 or 0 for true or false whether to plot the timecourse
% cell : which cell is the STA being computed for, necessary to figure out
% which part of the movie to use
% num_cells : total number of cells targeted with this movie so that it
% knows how many sections to parse the movie into


%% Results
% sta : Computed in matlab not java, binning problems
% timecourse: From sta sig stixels
% significant stixels : standard parameters

%% Author
% Colleen Rhoades (rhoades@stanford.edu)
% April 7, 2015


function [sta] = compute_jitter_sta_20160217(datarun, mdf_file, num_frames, spikes, jitter_x, jitter_y,  stixel_size, num_colors,save_path)
%% This function computes the STA without relying on STAs from vision. The binning is slightly different from Vision.
%     mglBltTexture(frametex, [stimulus.x_start+jitterX, stimulus.y_start+jitterY, stimulus.span_width, stimulus.span_height], -1, -1);
dbstop if error
% jitter_x = zeros(size(jitter_x));
% jitter_y = zeros(size(jitter_y));
%% ---------------------------------- Process movie ------------------------------
%onsets of the stimulus presentation
datarun.triggers = [datarun.triggers; [datarun.triggers(end) + mean(diff(datarun.triggers)):mean(diff(datarun.triggers)):datarun.triggers(end) + 300*mean(diff(datarun.triggers))]'];
%
datarun.triggers= datarun.triggers(1:300);
triggers = datarun.triggers;
% [mov,height,width,duration,refresh] = get_movie_ath(mdf_file,...
%     triggers, 1,2);

bw = 0;
[inputs, refresh, duration] = get_wn_movie_ath(datarun, mdf_file, bw);
% test = reshape(1:800*85000, 800, 85000);

% test = [repmat([1;2;3],800,1), repmat([4;5;6],800,1),repmat([7;8;9],800,1),repmat([10;11;12],800,1)]
% inputs = test;
% real_frame(:,1) = inputs(1:(25));


%% account for dropped frames
%% data026
%% 1857 is x coordinate of sharp peak in diff(triggers) graph
% triggers = [triggers(1:1857); triggers(1858:end) - mean(diff(triggers(1:1857)))];
% jitter_x_new = [jitter_x(1:92851); jitter_x(92900:end)]; %1857*100/2 (refresh is 2)+1 : 1857*100/2 (refresh is 2)+50
% jitter_y_new = [jitter_y(1:92851); jitter_y(92900:end)];
% inputs = [inputs(:,1:92851), inputs(:,92900:end)];
for i = 1:size(spikes,2)
    ind = find(spikes{i} > 92851/60 & spikes{i} <= 92900/60);
    if ~isempty(ind)
        spikes{i}(ind(end)+1:end) = spikes{i}(ind(end)+1:end) - (92900/60-92851/60);
        spikes{i} = [spikes{i}(1:ind(1)-1); spikes{i}(ind(end)+1:end)];
    end
    
    
end
%% data024
% triggers = [triggers(1:392); triggers(393:end) - mean(diff(triggers(1:392)))];
 jitter_x_new = [jitter_x]; %392*100/2 (refresh is 2)+1 : 392*100/2 (refresh is 2)+50
 jitter_y_new = [jitter_y];
% inputs = [inputs(:,1:19601), inputs(:,19650:end)];
% for i = 1:size(spikes,2)
%     ind = find(spikes{i} > 19601/60 & spikes{i} <= 19650/60);
%     spikes{i}(ind(2)+1:end) = spikes{i}(ind(2)+1:end) - (19650/60-19601/60);
%     spikes{i} = [spikes{i}(1:ind(1)-1); spikes{i}(ind(2)+1:end)];
%
% end




image_width = 40;
image_height = 20;

real_frame = zeros(image_width, image_height, 3);
real_frame(:,:,1,1) = reshape(inputs(1:3:image_width*image_height*3)',image_width, image_height);
real_frame(:,:,2,1) = reshape(inputs(2:3:image_width*image_height*3)',image_width, image_height);
real_frame(:,:,3,1) = reshape(inputs(3:3:image_width*image_height*3)',image_width, image_height);

pointer = image_width*image_height*3+1;
%     pointer = 2+25+2;
i =2;
while pointer+2*3+image_height*image_width*3-1<size(inputs,2)*size(inputs,1)
    temp = inputs(pointer+2*3:pointer+2*3+image_height*image_width*3-1);
    real_frame(:,:,1,i) = reshape(temp(1:3:end), image_width, image_height);
    real_frame(:,:,2,i) = reshape(temp(2:3:end), image_width, image_height);
    real_frame(:,:,3,i) = reshape(temp(3:3:end), image_width, image_height);
    
    pointer = pointer+2*3+image_height*image_width*3;
    i = i+1;
end

% a = reshape(real, 20, 40, size(real,2));

%inputs is 800x85000

% [mvi] = load_movie(mdf_file, triggers);
% mvi = squeeze(mvi(:,:,1,:));
% Compute the time each stimulus frame occurred
length_of_time = ceil(triggers(end))+1;
upsampled_num_frames = length_of_time*120;

upsample_factor = round(refresh/(100/12));
bt_triggers = triggers - [0;triggers(1:end-1)];
avg_bt_triggers = mean(bt_triggers(2:end));
triggers = [triggers; triggers(end) + avg_bt_triggers];
for j = 1:size(spikes,2)
    
    
    frames_needed{j} = zeros(3,(length(triggers)-2)*100+120);
    a = kron(1:upsampled_num_frames/upsample_factor, ones(1,upsample_factor));
    frames_needed{j}(1,:) = a(1:(length(triggers)-2)*100+120);
    
    %     frames_needed{j}(1,:) = kron(1:upsampled_num_frames/upsample_factor, ones(1,upsample_factor));
    
    
    for i= 1:length(triggers)-1
        spacing = linspace(triggers(i), triggers(i+1),101);
        frames_needed{j}(2, (i-1)*100+1:(i-1)*100+100)= spacing(1:end-1);
    end
    
    for i = 1:size(frames_needed{j},2)-1
        frames_needed{j}(3,i) = sum(spikes{j} >= frames_needed{j}(2,i) & spikes{j} < frames_needed{j}(2,i+1));
    end
end


%% Compute movie
movie = zeros(40*stixel_size, 20*stixel_size, 3,500);

sta = cell(size(spikes,2),1);
for i = 1:size(spikes,2)
    sta{i} =zeros(size(movie,1),size(movie,2),num_colors, num_frames);
end


start_points = [1:500:size(frames_needed{1},2) size(frames_needed{1},2)];

movie_exist = exist(['/Volumes/Lab/Users/crhoades/JitterMovie/2016-02-17-6/data026/', 'movie_block_1.mat']);
if movie_exist ==1000
    for j = 1:length(start_points)-1
        
        for m = 1:10000/200
            temp = load(['/Volumes/Lab/Users/crhoades/JitterMovie/2016-02-17-6/data026/', 'movie_block_', num2str(50*(j-1)+m)], 'current_movie');
            movie(:,:,:,200*(m-1)+1:200*(m-1)+200) = temp.current_movie;
        end
        
        %             current_movie = movie(:,:,:,200*(m-1)+1:200*(m-1)+200);
        %             save([save_path, 'movie_block_', num2str(50*(j-1)+m)], 'current_movie');
        for cel = 1:size(spikes,2)
            for i = 1:start_points(j+1)-1 - start_points(j)
                if mod(i,1000) == 0
                    fprintf('Phase: %d out of %d, %d out of %d \n', j, length(start_points)-1, i, start_points(j+1)-1 - start_points(j));
                end
                
                
                %         if start_points(j) -1+ i <= length(spikes_by_frame)
                if frames_needed{cel}(3,start_points(j)-1 + i) == 0
                else
                    if i <= num_frames
                    else
                        for t = 1:num_frames
                            subtract = num_frames - t +1;
                            sta{cel}(:,:, :,subtract) = sta{cel}(:,:, :, subtract) + movie(:,:,:,i-t) * frames_needed{cel}(3,start_points(j)-1 + i);
                        end
                    end
                end
            end
            
        end
        if ~exist(['/Volumes/Lab/Users/crhoades/Jitter/2016-02-17-6/data026cf_spikes'])
            mkdir(['/Volumes/Lab/Users/crhoades/Jitter/2016-02-17-6/data026cf_spikes']);
        end
        if mod(j,10) == 0
            
            save(['/Volumes/Lab/Users/crhoades/Jitter/2016-02-17-6/data026cf_spikes/temp'], 'sta', '-v7.3');
        end
    end
    
    
else
    height = 20;
    width = 40;
    for j = 1:length(start_points)-1
        for i = 1:start_points(j+1)-1 - start_points(j)
            
            if mod(i,1000) == 0
                fprintf('Phase: %d out of %d, %d out of %d \n', j, length(start_points)-1, i, start_points(j+1)-1 - start_points(j));
            end
            if start_points(j)-1 + i+num_frames < size(frames_needed{1},2) %&& start_points(j)_1 + i - num_frames>0
                %             if sum(frames_needed(3,(start_points(j) + i:start_points(j) + i-1+num_frames)))~=0
                %                 if i  <= (duration - 1)
                true_frame = zeros(width*stixel_size, height*stixel_size);
                F = real_frame(:,:,:,frames_needed{1}(1,start_points(j)-1 + i));
                %                     F = round(mvi.getFrame(frames_needed{cel}(1,start_points(j)-1 + i)).getBuffer);
                shaped_frame = F(:,:,1);
                sized_frame = imresize(double(shaped_frame), stixel_size, 'nearest');
                movie(:,:,1,i) = sized_frame;
                sized_frame = sized_frame((stixel_size/2+1):(end - stixel_size/2), (stixel_size/2+1):(end - stixel_size/2));
                position = [jitter_x_new(frames_needed{1}(1,start_points(j)-1 + i))+1+stixel_size/2, jitter_y_new(frames_needed{1}(1,start_points(j)-1 + i))+1+stixel_size/2];
                %         x and y might be reversed
                true_frame(position(1):(size(sized_frame,1)+position(1)-1), position(2):(size(sized_frame,2)+position(2)-1)) = sized_frame;
                movie(:,:,1,i) = true_frame;
                if num_colors == 3
                    shaped_frame = F(:,:,2);
                    sized_frame = imresize(double(shaped_frame), stixel_size, 'nearest');
                    sized_frame = sized_frame((stixel_size/2+1):(end - stixel_size/2), (stixel_size/2+1):(end - stixel_size/2));
                    % x and y might be reversed
                    true_frame(position(1):(size(sized_frame,1)+position(1)-1), position(2):(size(sized_frame,2)+position(2)-1)) = sized_frame;
                    movie(:,:,2,i) =true_frame;
                    
                    shaped_frame = F(:,:,3);
                    sized_frame = imresize(double(shaped_frame), stixel_size, 'nearest');
                    sized_frame = sized_frame((stixel_size/2+1):(end - stixel_size/2), (stixel_size/2+1):(end - stixel_size/2));
                    % x and y might be reversed
                    true_frame(position(1):(size(sized_frame,1)+position(1)-1), position(2):(size(sized_frame,2)+position(2)-1)) = sized_frame;
                    movie(:,:,3,i) = true_frame;
                end
                %                 else
                %                     continue
                %                 end
                %             else
                %                 counter = counter +1;
                %             end
                %
                
            end
            for cel = 1:size(spikes,2)
                
                %         if start_points(j) -1+ i <= length(spikes_by_frame)
                if frames_needed{cel}(3,start_points(j)-1 + i) == 0
                else
                    if i <= num_frames
                    else
                        for t = 1:num_frames
                            subtract = num_frames - t +1;
                            sta{cel}(:,:, :,subtract) = sta{cel}(:,:, :, subtract) + movie(:,:,:,i-t) * frames_needed{cel}(3,start_points(j)-1 + i);
                        end
                    end
                end
                %         end
                
            end
            
        end
        
        % save('/Volumes/Lab/Users/crhoades/Colleen/matlab/private/colleen/New Cell Types/Jitter/temp_sta', 'sta');
        
%         for m = 1:10000/200
%             current_movie = movie(:,:,:,200*(m-1)+1:200*(m-1)+200);
%             if ~exist(save_path)
%                 mkdir(save_path)
%             end
%             save([save_path, 'movie_block_', num2str(50*(j-1)+m)], 'current_movie');
%         end
        
        if ~exist(['/Volumes/Lab/Users/crhoades/Jitter/2016-02-17-6/data026cf_spikes'])
            mkdir(['/Volumes/Lab/Users/crhoades/Jitter/2016-02-17-6/data026cf_spikes']);
        end
        if mod(j,5) == 0
            
            save(['/Volumes/Lab/Users/crhoades/Jitter/2016-02-17-6/data026cf_spikes/temp'], 'sta', '-v7.3');
        end
    end
end

%     save(['/Volumes/Lab/Users/crhoades/JitterMovie/2008-04-22-5/data004/chuck_', num2str(counter)],'movie', '-v7.3');



%% ---------------------- Use spike times to form STA --------------------------







figure;
for i = 1:num_frames
    if size(sta{1},3) == 3
        imagesc(sta{1}(10:end-10,10:end-10,2,i))
    else
        imagesc(sta{1}(:,:,1,i))
    end
    pause(0.25)
end
%
%
% axis equal
% % sta = permute(sta,[1,2,4,3]);
% sig_stixels = significant_stixels(sta);
% rf = rf_from_sta(sta, 'sig_stixels', sig_stixels);
% figure;
% imagesc(norm_image(rf));
% axis equal







