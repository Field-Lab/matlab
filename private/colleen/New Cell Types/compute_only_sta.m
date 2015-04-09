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


function [sta, timecourse, sig_stixels] = compute_only_sta(datarun, mdf_file, num_frames, spikes, plotting, cell, num_cells)
%% This function computes the STA without relying on STAs from vision. The binning is slightly different from Vision.


%% ---------------------------------- Process movie ------------------------------
triggers=datarun.triggers; %onsets of the stimulus presentation

[mov,height,width,duration,refresh] = get_movie_ath(mdf_file,...
    triggers, 1,2);

[mvi] = load_movie(mdf_file, triggers);

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
frame_times = frame_times*1000; % in ms


height = height/num_cells; % assume movie file has the STA stacked vertically
% initialize STA
sta=zeros(height,width,3, num_frames); %height, width, frames back
tic
icnt=0;
spikes = spikes(1:10000); % only use for testing parasols

%% ---------------------- Use spike times to form STA --------------------------
for i=spikes'
    % ignore spikes without num_frames preceding it
    start=find(frame_times>i,1)-num_frames;
    if(start>000) % don't use the spikes that don't have num_frames before it
        icnt=icnt+1;
        if mod(icnt,1000) == 0
            fprintf('%d out of %d \n', icnt, length(spikes)')
        end
        
        for j=1:num_frames
            F = round(mvi.getFrame(start+j).getBuffer);
            pixels = length(F)/num_cells;
            sta(:,:,1, j) = sta(:,:,1,j) + round(reshape(F(cell*pixels-pixels+1:3:cell*pixels),width,height)'-0.5); % store the three color channels
            sta(:,:,2, j) = sta(:,:,2,j) + round(reshape(F(cell*pixels-pixels+2:3:cell*pixels),width,height)'-0.5);
            sta(:,:,3, j) = sta(:,:,3,j) + round(reshape(F(cell*pixels-pixels+3:3:cell*pixels),width,height)'-0.5);
            %             sta(:,:,1, j) = sta(:,:,1,j) + round(reshape(F(1:3:end),width,height)'-0.5); % store the three color channels
            %             sta(:,:,2, j) = sta(:,:,2,j) + round(reshape(F(2:3:end),width,height)'-0.5);
            %             sta(:,:,3, j) = sta(:,:,3,j) + round(reshape(F(3:3:end),width,height)'-0.5);
        end
    end
end
sta=sta/icnt;

% find best STA frame
[junk,start_index] = max(sum(reshape(sta.^2,[],size(sta,4)),1));

% normalize STA color
% sta = norm_image(sta);

%% ------------------ Compute the timecourse --------------------------
[sig_stixels] = significant_stixels(sta);
[timecourse, params] = time_course_from_sta(sta, sig_stixels);


if plotting == 1    
    h = plot_time_course_(timecourse, 'colors', ['rgb']', 'foa', 0)
    title('TimeCourse from MATLAB STA')
    % Show the best frame of the STA
    figure
    imagesc(squeeze(norm_image(sta(:,:,:, start_index))));
    title(['STA Frame: ' num2str(start_index)]);
end
















