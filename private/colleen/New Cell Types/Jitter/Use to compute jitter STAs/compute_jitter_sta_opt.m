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


function [sta] = compute_jitter_sta_opt(datarun, mdf_file, num_frames, spikes, jitter_x, jitter_y,  stixel_size, num_colors,dataparam)
%% This function computes the STA without relying on STAs from vision. The binning is slightly different from Vision.
%     mglBltTexture(frametex, [stimulus.x_start+jitterX, stimulus.y_start+jitterY, stimulus.span_width, stimulus.span_height], -1, -1);
dbstop if error

%% ---------------------------------- Process movie ------------------------------

triggers = datarun.triggers;

if num_colors == 3
    bw = 0;
else
    bw= 1;
end

[inputs, ~, ~] = get_wn_movie_ath(datarun, mdf_file, bw);


%% account for dropped frames

difference = diff(triggers);
mean_val = mean(difference);
high_point_trigger = find(difference>(mean_val*1.2)==1); 
high_point_trigger = sort(high_point_trigger,'descend');
for i = 1:length(high_point_trigger)
    triggers = [datarun.triggers(1:high_point_trigger(i)); datarun.triggers((high_point_trigger(i)+1):end) - mean(diff(datarun.triggers(1:high_point_trigger(i))))];
    %1496*100/refresh (refresh is 2)+1 : 392*100/2 (refresh is 2)+50
    
    % inputs = [inputs(:,1:149601), inputs(:,(149601 + (100/dataparam.interval)-1):end)];
    for j = 1:size(spikes,2)
        ind = find(spikes{j} > datarun.triggers(high_point_trigger(i)) & spikes{j} <= (datarun.triggers(high_point_trigger(i)+1)+num_frames*mean(diff(datarun.triggers(1:high_point_trigger(i))))/(100/dataparam.interval))); % rate = 120 if refresh =1 or 60 if refresh = 2
        if ~isempty(ind)
            spikes{j}(ind(end)+1:end) = spikes{j}(ind(end)+1:end) - mean(diff(datarun.triggers(1:high_point_trigger(i))));
            spikes{j} = [spikes{j}(1:ind(1)-1); spikes{j}(ind(end)+1:end)];
        end
    end
    
end



image_width = dataparam.x_dim/stixel_size;
image_height = dataparam.y_dim/stixel_size;

real_frame = zeros(image_width, image_height, num_colors, size(inputs,2), 'int16');
if num_colors == 3
    real_frame(:,:,1,1) = reshape(inputs(1:3:image_width*image_height*3)',image_width, image_height);
    real_frame(:,:,2,1) = reshape(inputs(2:3:image_width*image_height*3)',image_width, image_height);
    real_frame(:,:,3,1) = reshape(inputs(3:3:image_width*image_height*3)',image_width, image_height);
else
    real_frame(:,:,1,1) = reshape(inputs(1:image_width*image_height)',image_width, image_height);
    
end


pointer = image_width*image_height*num_colors+1;
i =2;

if strcmp(dataparam.date(1:10), '2016-02-17') % bug in jitter code on this date
    %%%% REMOVE THE +6 WHEN NOT 2016-02-17

    while pointer+2*3+image_height*image_width*3-1<size(inputs,2)*size(inputs,1)
        temp = inputs(pointer+2*3:pointer+2*3+image_height*image_width*3-1);
        real_frame(:,:,1,i) = reshape(temp(1:3:end), image_width, image_height);
        real_frame(:,:,2,i) = reshape(temp(2:3:end), image_width, image_height);
        real_frame(:,:,3,i) = reshape(temp(3:3:end), image_width, image_height);

        pointer = pointer+2*3+image_height*image_width*3;
        i = i+1;
    end
else
    
    while pointer+image_height*image_width*num_colors-1<size(inputs,2)*size(inputs,1)
    temp = inputs(pointer:pointer+image_height*image_width*num_colors-1);
    if num_colors == 3
        real_frame(:,:,1,i) = reshape(temp(1:3:end), image_width, image_height);
        real_frame(:,:,2,i) = reshape(temp(2:3:end), image_width, image_height);
        real_frame(:,:,3,i) = reshape(temp(3:3:end), image_width, image_height);
    else
        real_frame(:,:,1,i) = reshape(temp, image_width, image_height);
    end
    
    pointer = pointer+image_height*image_width*num_colors;
    i = i+1;
    end


end



upsampled_num_frames = length(triggers)*100/dataparam.interval;

upsample_factor = dataparam.interval; % should be interval

frames_needed = kron(1:ceil(upsampled_num_frames), ones(1,upsample_factor));

frame_spacing = zeros(1, size(frames_needed,2));
for i= 1:length(triggers)-1
    spacing = linspace(triggers(i), triggers(i+1),100 +1);
    frame_spacing(1, (i-1)*100+1:(i-1)*100+100)= spacing(1:end-1); %% assume triggers every 100 frames
end

binned_spikes = zeros(size(spikes,2), size(frames_needed,2)-1, 'int16');
for j = 1:size(spikes,2)
    for i = 1:size(frames_needed,2)-1
        binned_spikes(j,i) = sum(spikes{j} >= frame_spacing(1,i) & spikes{j} < frame_spacing(1,i+1));
    end
end

%% Compute movie

fprintf('Progress: %s%s\n', dataparam.date, dataparam.concatname)
fprintf([repmat('-', 1,dataparam.num_of_interval+1), '\n'])
fprintf('*')
sta = zeros(image_width*stixel_size, image_height*stixel_size,num_colors, num_frames, size(spikes,2), 'int16');
sum_binned_spikes = sum(binned_spikes,1 );
try
    height =image_height;
    width = image_width;
    
    for i = 1:size(frames_needed,2)
        if mod(i,floor(size(frames_needed,2)/dataparam.num_of_interval)) == 0
            fprintf('*');
        end
        
        
        
        % don't compute the frame if you don't need it
        
        if sum(sum_binned_spikes(:,i + 1:i+num_frames)) ~= 0
          
            movie = zeros(image_width*stixel_size, image_height*stixel_size, num_colors, 'int16');
            true_frame = zeros(width*stixel_size, height*stixel_size, 'int16');
            F = real_frame(:,:,:,frames_needed(1,i));
            shaped_frame = F(:,:,1);
            %                 sized_frame = imresize(shaped_frame, stixel_size, 'nearest');
            
            scale = [stixel_size stixel_size]; % The resolution scale factors: [rows columns]
            oldSize = size(shaped_frame); % Get the size of your image
            newSize = scale.*oldSize;  % Compute the new image size
            
            % Compute an upsampled set of indices:
            
            rowIndex = min(round(((1:newSize(1))-0.5)./scale(1)+0.5),oldSize(1));
            colIndex = min(round(((1:newSize(2))-0.5)./scale(2)+0.5),oldSize(2));
            % Index old image to get new image:
            sized_frame = shaped_frame(rowIndex,colIndex);
            
            
            
            sized_frame = sized_frame((stixel_size/2+1):(end - stixel_size/2), (stixel_size/2+1):(end - stixel_size/2));
            position = [jitter_x(frames_needed(1,i))+1+stixel_size/2, jitter_y(frames_needed(1,i))+1+stixel_size/2];
            true_frame(position(1):(size(sized_frame,1)+position(1)-1), position(2):(size(sized_frame,2)+position(2)-1)) = sized_frame;
            movie(:,:,1) = true_frame;
            if num_colors == 3
                shaped_frame = F(:,:,2);
                sized_frame = shaped_frame(rowIndex,colIndex);
                sized_frame = sized_frame((stixel_size/2+1):(end - stixel_size/2), (stixel_size/2+1):(end - stixel_size/2));
                true_frame(position(1):(size(sized_frame,1)+position(1)-1), position(2):(size(sized_frame,2)+position(2)-1)) = sized_frame;
                movie(:,:,2) =true_frame;
                
                shaped_frame = F(:,:,3);
                sized_frame = shaped_frame(rowIndex,colIndex);
                sized_frame = sized_frame((stixel_size/2+1):(end - stixel_size/2), (stixel_size/2+1):(end - stixel_size/2));
                true_frame(position(1):(size(sized_frame,1)+position(1)-1), position(2):(size(sized_frame,2)+position(2)-1)) = sized_frame;
                movie(:,:,3) = true_frame;
            end
            
            
            
            
            
            
            
            for t = 1:num_frames

                [x_ind] = find(binned_spikes(:,i+t)>0);
                if ~isempty(x_ind)                        
                    subtract = num_frames -t+1;

                    if subtract > 0
                        for cel = 1:length(x_ind)
                            % maximum values in sta are -32768 and
                            % 32768, which could be problem for cells
                            % that spike A LOT
                            sta(:,:,:,subtract, x_ind(cel)) = sta(:,:,:,subtract, x_ind(cel))  + movie * binned_spikes(x_ind(cel),i+t);
                        end
                    end
                    
                end
            end
        end
        
    end
catch
    disp('out of frames')
end
fprintf('\n');

sta = double(sta);
for i = 1:size(binned_spikes,1)
    sta(:,:,:,:,i) = sta(:,:,:,:,i)./sum(double(binned_spikes(i,:)));
end









