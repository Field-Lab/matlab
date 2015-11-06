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
% frame_times = frame_times*1000; % in ms


% height = height/num_cells; % assume movie file has the STA stacked vertically
% initialize STA
sta=zeros(height*stixel_size,width*stixel_size,num_frames, 3); %height, width, frames back
% sta=zeros(height,width,num_frames, 3); %height, width, frames back

tic
icnt=0;
% spikes = spikes(1:10000); % only use for testing parasols
% spikes = spikes(1:4199);
if length(spikes) > 5000
    spikes = spikes(1:5000);
end

%% Compute movie
movie = zeros(height*stixel_size, width*stixel_size, duration, 3);
% counter = 1;
% for iter= 1:600:duration-1
    for i = 1:duration%600
        if i  <= (duration - 1)
            true_frame = zeros(height*stixel_size, width*stixel_size);

            F = round(mvi.getFrame(i).getBuffer);
            shaped_frame = round(reshape(F(1:3:end),width,height)'-0.5);
            sized_frame = imresize(double(shaped_frame), stixel_size, 'nearest');
            %     movie(:,:,i,1) = sized_frame;
            sized_frame = sized_frame((stixel_size/2+1):(end - stixel_size/2), (stixel_size/2+1):(end - stixel_size/2));
            position = [jitter_x(i)+1+stixel_size/2, jitter_y(i)+1+stixel_size/2];
            % x and y might be reversed
            true_frame(position(1):(size(sized_frame,1)+position(1)-1), position(2):(size(sized_frame,2)+position(2)-1)) = sized_frame;
            movie(:,:,i,1) = true_frame;
            if num_colors == 3
                shaped_frame = round(reshape(F(2:3:end),width,height)'-0.5);
                sized_frame = imresize(double(shaped_frame), stixel_size, 'nearest');
                sized_frame = sized_frame((stixel_size/2+1):(end - stixel_size/2), (stixel_size/2+1):(end - stixel_size/2));
                position = [jitter_x(i)+1+stixel_size/2, jitter_y(i)+1+stixel_size/2];
                % x and y might be reversed
                true_frame(position(1):(size(sized_frame,1)+position(1)-1), position(2):(size(sized_frame,2)+position(2)-1)) = sized_frame;
                movie(:,:,i,2) = true_frame;

                shaped_frame = round(reshape(F(3:3:end),width,height)'-0.5);
                sized_frame = imresize(double(shaped_frame), stixel_size, 'nearest');
                sized_frame = sized_frame((stixel_size/2+1):(end - stixel_size/2), (stixel_size/2+1):(end - stixel_size/2));
                position = [jitter_x(i)+1+stixel_size/2, jitter_y(i)+1+stixel_size/2];
                % x and y might be reversed
                true_frame(position(1):(size(sized_frame,1)+position(1)-1), position(2):(size(sized_frame,2)+position(2)-1)) = sized_frame;
                movie(:,:,i,3) = true_frame;
            else
                 movie(:,:,i,2) = true_frame;
                 movie(:,:,i,3) = true_frame;
            end
        else
            continue
        end
        
    end
    
%     save(['/Volumes/Lab/Users/crhoades/JitterMovie/2008-04-22-5/data004/chuck_', num2str(counter)],'movie', '-v7.3');



%% ---------------------- Use spike times to form STA --------------------------
   
    
    spikes_by_frame = nan(length(frame_times)-1,1);
    for i = 1:length(frame_times)-1
        spikes_by_frame(i) = sum(spikes >= frame_times(i) & spikes < frame_times(i+1));
    end
    frameR = nan(height*width*stixel_size*stixel_size, length(frame_times)-1);
%     frameG = nan(height*width*stixel_size*stixel_size, length(frame_times)-1);
%     frameB= nan(height*width*stixel_size*stixel_size,length(frame_times)-1);
    
    for i = 1:length(frame_times)-1
        frameR(:, i) = reshape(squeeze(movie(:,:,i,1)), [size(squeeze(movie(:,:,i,1)), 1)*size(squeeze(movie(:,:,i,1)), 2), 1]);
%         frameG(:,i) = reshape(squeeze(movie(:,:,i,2)), [size(squeeze(movie(:,:,i,2)), 1)*size(squeeze(movie(:,:,i,2)), 2), 1]);
%         frameB(:, i) = reshape(squeeze(movie(:,:,i,3)), [size(squeeze(movie(:,:,i,3)), 1)*size(squeeze(movie(:,:,i,3)), 2), 1]);
        
    end
    
%     sta= nan(height,width, 3,num_frames);
    for i= 1:num_frames
        subtract= num_frames -i +1;
        one_frameR = frameR(:,1:end-subtract+1)*spikes_by_frame(subtract:end);
%         one_frameG = frameG(:,1:end-subtract+1)*spikes_by_frame(subtract:end);
%         one_frameB = frameB(:,1:end-subtract+1)*spikes_by_frame(subtract:end);
        
        sta(:,:,i,1) = reshape(one_frameR, width, height)';
%         sta{cells}(:,:,2,i) = reshape(one_frameG, width, height)';
%         sta{cells}(:,:,3,i) = reshape(one_frameB, width, height)';
        
    end
    
    sta{cells} = sta{cells}/length(spikes);
    
% for i=spikes'
%     % ignore spikes without num_frames preceding it
%     start=find(frame_times>i,1)-num_frames;
%     if(start>000) % don't use the spikes that don't have num_frames before it
%         icnt=icnt+1;
%         if mod(icnt,1000) == 0
%             fprintf('%d out of %d \n', icnt, length(spikes)')
%         end
%         
%         for j=1:num_frames
%             try
%                 
%                 F = movie(:,:,start+j, :);
%                 sta(:,:,j,1) = sta(:,:,j,1) + F(:,:,1,1); % store the three color channels
%                 
%                 if num_colors == 3
%                     sta(:,:,j,2) = sta(:,:,j,2) +  F(:,:,1,2); % store the three color channels
%                     
%                      sta(:,:,j,3) = sta(:,:,j,3) +  F(:,:,1,3); % store the three color channels
%                     
%                 end
%                  catch
%                 a = 1;
%             end
%             
%         end
%     end
% end
% if num_colors ==1
%     sta = sta(:,:,:,1);
% end

% % need to rearrange the sta to be in the map order
% sta_new = zeros(size(map,1)/stixel_size,size(map,2)/stixel_size,num_colors,num_frames);
% map_x = map(:,1:stixel_size:end);
% map_stix = map_x(1:stixel_size:end,:);
% % sta = permute(sta, [2,1,3,4]); % so that single indexing goes vertically
% for i = 1:height*width
%     [x,y] = find(map_stix == i);
%     down_index = floor((i-1)/width)
%     across_index = mod(i-1,width)
%     sta_new(x,y, :,:) = sta(down_index+1,across_index+1,:,:);
% end
sta=sta/icnt;
% figure;
% for i =27%:30
% imagesc(sta(:,:,i,1))
% axis equal
% pause(0.5)
% end

% for i = 24:30
%     if num_colors == 3
%         figure; imagesc(sta_new(:,:,2,i))
%     else
%         figure; imagesc(sta_new(:,:,1,i))
%     end
%     
% end


% image(offset_y+1: offset_y+height, offset_x+1:offset_x+width, :,:) = sta;

% sta = image;
% find best STA frame
% [junk,start_index] = max(sum(reshape(sta_new.^2,[],size(sta_new,4)),1));
% 
% % normalize STA color
% % sta = norm_image(sta);
% close al
% %% ------------------ Compute the timecourse --------------------------
% [sig_stixels] = significant_stixels(sta_new);
% 
% [timecourse, params] = time_course_from_sta(sta_new, [18,24])%sig_stixels);
% 
% 
% if plotting == 1
%     h = plot_time_course_(timecourse, 'colors', ['rgb']', 'foa', 0)
%     title('TimeCourse from MATLAB STA')
%     % Show the best frame of the STA
%     figure
%     imagesc(squeeze(norm_image(sta_new(:,:,:, start_index))));
%     title(['STA Frame: ' num2str(start_index)]);
% end
















