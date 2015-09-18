%% --------------- Compute STA of focal white noise in MATLAB (not Java) -------------------
%% Function: Compute STA in matlab directly

function [sta, timecourse, sig_stixels] = compute_matlab_sta_fast(datarun, dataparam, vision_id)
%% This function computes the STA without relying on STAs from vision. The binning is slightly different from Vision.
mdf_file = dataparam.mdf_file;
num_frames = dataparam.num_frames;
cell_index = get_cell_indices(datarun, vision_id);
to_save = dataparam.to_save;
save_location = ['/Users/colleen/Desktop/MATLAB_STAs/', dataparam.date, '/', dataparam.concatname, '/'];
if exist('dataparam.scale')
    scale = dataparam.scale;
else
    scale = 10;
end

if ~exist(save_location)
    mkdir(save_location);
end


%% ---------------------------------- Process movie ------------------------------
triggers=datarun.triggers; %onsets of the stimulus presentation

[mov,height,width,duration,refresh] = get_movie_ath(mdf_file,...
    triggers, 1,2);

[mvi] = load_movie(mdf_file, triggers);

bt_triggers = triggers - [0;triggers(1:end-1)];
avg_bt_triggers = mean(bt_triggers(2:end));
frames_per_trigger = round(avg_bt_triggers*1000/refresh);
last_trigger_time = ceil(triggers(end));
frame_times = zeros(ceil(last_trigger_time/refresh*1000),1);
for i = 1: length(triggers)
    temp = linspace(triggers(i),triggers(i)+ (frames_per_trigger-1)*refresh/1000,frames_per_trigger)';
    frame_times(i*frames_per_trigger-(frames_per_trigger-1):i*frames_per_trigger) = temp;
end
sta = cell(length(cell_index),1);
for cells = 1:length(cell_index)
    spikes = datarun.spikes{cell_index(cells)};
    
    
    spikes_by_frame = nan(length(frame_times)-1,1);
    for i = 1:length(frame_times)-1
        spikes_by_frame(i) = sum(spikes >= frame_times(i) & spikes < frame_times(i+1));
    end
    frameR = nan(height*width, length(frame_times)-1);
    frameG = nan(height*width, length(frame_times)-1);
    frameB= nan(height*width,length(frame_times)-1);
    
    for i = 1:length(frame_times)-1
        F = round(mvi.getFrame(i).getBuffer);
        frameR(:, i) = round(F(1:3:end)-0.5);
        frameG(:,i) = round(F(2:3:end)-0.5);
        frameB(:, i) = round(F(3:3:end)-0.5);
        
    end
    
    sta{cells} = nan(height,width, 3,num_frames);
    for i= 1:num_frames
        subtract= num_frames -i +1;
        one_frameR = frameR(:,1:end-subtract+1)*spikes_by_frame(subtract:end);
        one_frameG = frameG(:,1:end-subtract+1)*spikes_by_frame(subtract:end);
        one_frameB = frameB(:,1:end-subtract+1)*spikes_by_frame(subtract:end);
        
        sta{cells}(:,:,1,i) = reshape(one_frameR, width, height)';
        sta{cells}(:,:,2,i) = reshape(one_frameG, width, height)';
        sta{cells}(:,:,3,i) = reshape(one_frameB, width, height)';
        
    end
    
    sta{cells} = sta{cells}/length(spikes);
    upsample_sta{cells} = imresize(sta{cells}, scale,'nearest');
    
    figure;
    subplot(2,1,1)
    imagesc(norm_image(squeeze(upsample_sta{cells}(:,:,:,26))));
    hold on
    axis equal
    axis off
    [sig_stixels{cells}] = significant_stixels(upsample_sta{cells});
    
    boundary = bwboundaries(full(sig_stixels{cells}));
    plot(boundary{1}(:,2)+0.5, boundary{1}(:,1)+0.5, 'r','LineWidth',0.5);
    
    [timecourse{cells}, params] = time_course_from_sta(upsample_sta{cells},sig_stixels{cells});
    subplot(2,1,2)
    set(gca,'ColorOrder',[1 0 0;0 1 0;0 0 1])
    
    hold all
    plot(timecourse{cells})
    if to_save
        hgexport(gcf, [save_location, 'Cell_',num2str(vision_id(cells))])
        data.sta = sta{cells}
        [sig_stixels_real{cells}] = significant_stixels(sta{cells});
        [timecourse_real{cells}, params] = time_course_from_sta(sta{cells},sig_stixels_real{cells});
        
        data.sig_stixels = sig_stixels_real{cells};
        data.timecourse = timecourse_real{cells};
        
        save([save_location, 'Cell_',num2str(vision_id(cells))], 'data')
    end
    
end















