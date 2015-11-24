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


function [sta] = compute_jitter_sta(datarun, mdf_file, num_frames, spikes, jitter_x, jitter_y,  stixel_size, num_colors,save_path)
%% This function computes the STA without relying on STAs from vision. The binning is slightly different from Vision.
%     mglBltTexture(frametex, [stimulus.x_start+jitterX, stimulus.y_start+jitterY, stimulus.span_width, stimulus.span_height], -1, -1);


%% ---------------------------------- Process movie ------------------------------
triggers=datarun.triggers; %onsets of the stimulus presentation

[mov,height,width,duration,refresh] = get_movie_ath(mdf_file,...
    triggers, 1,2);

[mvi] = load_movie(mdf_file, triggers);
% mvi = squeeze(mvi(:,:,1,:));
% Compute the time each stimulus frame occurred
length_of_time = ceil(triggers(end))+1;
upsampled_num_frames = length_of_time*120;

upsample_factor = round(refresh/(100/12));
frames_needed = zeros(3,(length(triggers)-2)*100+100+120);
frames_needed(1,:) = kron(1:upsampled_num_frames/upsample_factor, ones(1,upsample_factor));
bt_triggers = triggers - [0;triggers(1:end-1)];
avg_bt_triggers = mean(bt_triggers(2:end));
triggers = [triggers; triggers(end) + avg_bt_triggers];

for i= 1:length(triggers)-1
    spacing = linspace(triggers(i), triggers(i+1),101);
    frames_needed(2, (i-1)*100+1:(i-1)*100+100)= spacing(1:end-1);
end

for i = 1:size(frames_needed,2)-1
    frames_needed(3,i) = sum(spikes >= frames_needed(2,i) & spikes < frames_needed(2,i+1));
end


%% Compute movie
movie = zeros(height*stixel_size, width*stixel_size, num_colors,duration);

sta =zeros(size(movie,1),size(movie,2),num_colors, num_frames);


start_points = [1:10000:size(frames_needed,2) size(frames_needed,2)];

movie_exist = exist([save_path, 'movie_block_1.mat']);
if movie_exist ~= 0
    for j = 1:length(start_points)-1
        
        for m = 1:10000/200
            temp = load([save_path, 'movie_block_', num2str(50*(j-1)+m)], 'current_movie');
            movie(:,:,:,200*(m-1)+1:200*(m-1)+200) = temp.current_movie;
        end
        
        %             current_movie = movie(:,:,:,200*(m-1)+1:200*(m-1)+200);
        %             save([save_path, 'movie_block_', num2str(50*(j-1)+m)], 'current_movie');
        
        for i = 1:start_points(j+1)-1 - start_points(j)
            if mod(i,1000) == 0
                fprintf('Phase: %d out of %d, %d out of %d \n', j, length(start_points)-1, i, start_points(j+1)-1 - start_points(j));
            end
            
            %         if start_points(j) -1+ i <= length(spikes_by_frame)
            if frames_needed(3,start_points(j)-1 + i) == 0
            else
                if i <= num_frames
                else
                    for t = 1:num_frames
                        subtract = num_frames - t +1;
                        sta(:,:, :,subtract) = sta(:,:, :, subtract) + movie(:,:,:,i-t) * frames_needed(3,start_points(j)-1 + i);
                    end
                end
            end
            %         end
            
        end
    end
    
    
else
    
    for j = 1:length(start_points)-1
        for i = 1:start_points(j+1)-1 - start_points(j)
            if mod(i,1000) == 0
                fprintf('Phase: %d out of %d, %d out of %d \n', j, length(start_points)-1, i, start_points(j+1)-1 - start_points(j));
            end
            if start_points(j)-1 + i+num_frames < size(frames_needed,2) %&& start_points(j)_1 + i - num_frames>0
                %             if sum(frames_needed(3,(start_points(j) + i:start_points(j) + i-1+num_frames)))~=0
                %                 if i  <= (duration - 1)
                true_frame = zeros(height*stixel_size, width*stixel_size);
                F = round(mvi.getFrame(frames_needed(1,start_points(j)-1 + i)).getBuffer);
                shaped_frame = round(reshape(F(1:3:end),width,height)'-0.5);
                sized_frame = imresize(double(shaped_frame), stixel_size, 'nearest');
                movie(:,:,1,i) = sized_frame;
                sized_frame = sized_frame((stixel_size/2+1):(end - stixel_size/2), (stixel_size/2+1):(end - stixel_size/2));
                position = [jitter_x(frames_needed(1,start_points(j)-1 + i))+1+stixel_size/2, jitter_y(frames_needed(1,start_points(j)-1 + i))+1+stixel_size/2];
                %         x and y might be reversed
                true_frame(position(1):(size(sized_frame,1)+position(1)-1), position(2):(size(sized_frame,2)+position(2)-1)) = sized_frame;
                movie(:,:,1,i) = true_frame;
                if num_colors == 3
                    shaped_frame = round(reshape(F(2:3:end),width,height)'-0.5);
                    sized_frame = imresize(double(shaped_frame), stixel_size, 'nearest');
                    sized_frame = sized_frame((stixel_size/2+1):(end - stixel_size/2), (stixel_size/2+1):(end - stixel_size/2));
                    % x and y might be reversed
                    true_frame(position(1):(size(sized_frame,1)+position(1)-1), position(2):(size(sized_frame,2)+position(2)-1)) = sized_frame;
                    movie(:,:,2,i) =true_frame;
                    
                    shaped_frame = round(reshape(F(3:3:end),width,height)'-0.5);
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
            %         if start_points(j) -1+ i <= length(spikes_by_frame)
            if frames_needed(3,start_points(j)-1 + i) == 0
            else
                if i <= num_frames
                else
                    for t = 1:num_frames
                        subtract = num_frames - t +1;
                        sta(:,:, :,subtract) = sta(:,:, :, subtract) + movie(:,:,:,i-t) * frames_needed(3,start_points(j)-1 + i);
                    end
                end
            end
            %         end
            
        end
        for m = 1:10000/200
            current_movie = movie(:,:,:,200*(m-1)+1:200*(m-1)+200);
            save([save_path, 'movie_block_', num2str(50*(j-1)+m)], 'current_movie');
        end
    end
end

%     save(['/Volumes/Lab/Users/crhoades/JitterMovie/2008-04-22-5/data004/chuck_', num2str(counter)],'movie', '-v7.3');



%% ---------------------- Use spike times to form STA --------------------------







figure;
for i = 1:num_frames
    if size(sta,3) == 3
        imagesc(sta(:,:,2,i))
    else
        imagesc(sta(:,:,1,i))
    end
    pause(0.25)
end


axis equal
% sta = permute(sta,[1,2,4,3]);
sig_stixels = significant_stixels(sta);
rf = rf_from_sta(sta, 'sig_stixels', sig_stixels);
figure;
imagesc(norm_image(rf));
axis equal







