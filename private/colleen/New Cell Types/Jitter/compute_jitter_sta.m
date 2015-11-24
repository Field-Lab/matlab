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


function [sta] = compute_jitter_sta(datarun, mdf_file, num_frames, spikes, jitter_x, jitter_y,  stixel_size, num_colors)
%% This function computes the STA without relying on STAs from vision. The binning is slightly different from Vision.
%     mglBltTexture(frametex, [stimulus.x_start+jitterX, stimulus.y_start+jitterY, stimulus.span_width, stimulus.span_height], -1, -1);


%% ---------------------------------- Process movie ------------------------------
triggers=datarun.triggers; %onsets of the stimulus presentation

[mov,height,width,duration,refresh] = get_movie_ath(mdf_file,...
    triggers, 1,2);

[mvi] = load_movie(mdf_file, triggers);
% mvi = squeeze(mvi(:,:,1,:));
% Compute the time each stimulus frame occurred

bt_triggers = triggers - [0;triggers(1:end-1)];
avg_bt_triggers = mean(bt_triggers(2:end));
frames_per_trigger = round(avg_bt_triggers*1000/refresh);
last_trigger_time = ceil(triggers(end));
frame_times = zeros(ceil(last_trigger_time/refresh*1000),1);
for i = 1: length(triggers)
    temp = linspace(triggers(i),triggers(i)+ (frames_per_trigger-1)*refresh/1000,frames_per_trigger)';
    frame_times(i*frames_per_trigger-(frames_per_trigger-1):i*frames_per_trigger) = temp;
end

% initialize STA
sta=zeros(height*stixel_size,width*stixel_size,num_frames, num_colors); %height, width, frames back

tic
icnt=0;

spikes_by_frame = nan(length(frame_times)-1,1);
for i = 1:length(frame_times)-1
    spikes_by_frame(i) = sum(spikes >= frame_times(i) & spikes < frame_times(i+1));
end


%% Compute movie
movie = zeros(height*stixel_size, width*stixel_size, num_colors,duration);
counter = 0;
sta =zeros(size(movie,1),size(movie,2),num_colors, num_frames);

segment = ceil(duration/10000);

start_points = floor(linspace(1,duration, segment+1))%1:segment:duration;

% end_points = linspace(segment,length(start_points)*segment,length(start_points));
% end_points(end) = duration;
for j = 1:length(start_points)-1
    for i = 1:start_points(j+1)-1 - start_points(j)
        if mod(i,1000) == 0
            fprintf('%d out of %d \n', i, duration);
        end
        if start_points(j)-1 + i < length(spikes_by_frame) && start_points(j)-1 + i - num_frames>0
            if sum(spikes_by_frame(start_points(j) + i-1-num_frames:start_points(j) + i-2))~=0
                if i  <= (duration - 1)
                    true_frame = zeros(height*stixel_size, width*stixel_size);
                    F = round(mvi.getFrame(start_points(j)-1 + i).getBuffer);
                    shaped_frame = round(reshape(F(1:3:end),width,height)'-0.5);
                    sized_frame = imresize(double(shaped_frame), stixel_size, 'nearest');
                    movie(:,:,1,i) = sized_frame;
                    sized_frame = sized_frame((stixel_size/2+1):(end - stixel_size/2), (stixel_size/2+1):(end - stixel_size/2));
                    position = [jitter_x(start_points(j)-1 + i)+1+stixel_size/2, jitter_y(start_points(j) -1 + i)+1+stixel_size/2];
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
                else
                    continue
                end
            else
                counter = counter +1;
            end
            
            
        end
        if start_points(j) -1+ i <= length(spikes_by_frame)
            if spikes_by_frame(start_points(j) -1+ i)  == 0
            else
                if i <= num_frames
                else
                    for t = 1:num_frames
                        subtract = num_frames - t +1;
                        sta(:,:, :,subtract) = sta(:,:, :, subtract) + movie(:,:,:,i-t) * spikes_by_frame(start_points(j)-1 + i);
                    end
                end
            end
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







