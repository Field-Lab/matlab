function [sta, timecourse] = compute_only_sta(datarun, mdf_file, num_frames, spikes, plotting)
%% This function computes the STA without relying on STAs from vision. The binning is slightly different from Vision.
%% Requires datarun for the trigger infomation, the mdf_file for the movie, num_frames, the spikes sequence, and a 0 or 1 for if you want plotting

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);
datarun=load_data(datarun,opt);

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




% initialize STA
sta=zeros(height,width,3, num_frames); %height, width, frames back
tic
icnt=0;

for i=spikes'
    % ignore spikes without num_frames preceding it
    start=find(frame_times>i,1)-num_frames;
    if(start>000)
        icnt=icnt+1;
        if mod(icnt,1000) == 0
           fprintf('%d out of %d \n', icnt, length(spikes)')
        end
        
        for j=1:num_frames
            F = round(mvi.getFrame(start+j).getBuffer);
            sta(:,:,1, j) = sta(:,:,1,j) + round(reshape(F(1:3:end),width,height)'-0.5); % store the three color channels
            sta(:,:,2, j) = sta(:,:,2,j) + round(reshape(F(2:3:end),width,height)'-0.5);
            sta(:,:,3, j) = sta(:,:,3,j) + round(reshape(F(3:3:end),width,height)'-0.5);
        end
    end
end
sta=sta/icnt;

% find best STA frame
[junk,start_index] = max(sum(reshape(sta.^2,[],size(sta,4)),1));

% normalize STA color
sta = norm_image(sta);

% Compute the timecourse
[sig_stixels] = significant_stixels(sta);
[timecourse, params] = time_course_from_sta(sta, sig_stixels);


if plotting == 1
    
    
    h = plot_time_course_(timecourse, 'colors', ['rgb']', 'foa', 0)
    % Show the best frame of the STA
    figure
    imagesc(squeeze(sta(:,:,:, start_index)));
    title(['STA Frame: ' num2str(start_index)]);
end
















