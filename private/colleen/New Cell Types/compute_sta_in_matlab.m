clear
%% Get timecourse of related cell

datarun.names.rrs_neurons_path='/Volumes/Analysis/2007-03-27-2/data003-gdf/data003.neurons';

mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-20-4-0.48-11111.xml';
num_frames = 14; % both have to be run with the name number of frames
target_cell = 7400;


opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);
datarun=load_data(datarun,opt);

triggers=datarun.triggers; %onsets of the stimulus presentation

[mov,height,width,duration,refresh] = get_movie_ath(mdf_file,...
    triggers, 1,2);

[mvi] = load_movie(mdf_file, triggers);

%frames per trigger

bt_triggers = triggers - [0;triggers(1:end-1)];
avg_bt_triggers = mean(bt_triggers(2:end));
frames_per_trigger = round(avg_bt_triggers*1000/refresh);

last_trigger_time = ceil(triggers(end));

frame_times = zeros(ceil(last_trigger_time/refresh*1000),1);

% frame_times(1:25:end) = triggers

temp = linspace(triggers(1),refresh/1000,frames_per_trigger);
for i = 1: length(triggers)
    temp = linspace(triggers(i),triggers(i)+ (frames_per_trigger-1)*refresh/1000,frames_per_trigger)';
    frame_times(i*frames_per_trigger-(frames_per_trigger-1):i*frames_per_trigger) = temp;
end
frame_times = frame_times*1000; % in ms

cellID=find(datarun.cell_ids==target_cell)

spikes=datarun.spikes{cellID};
% spikes are in s. convert to ms
spikes=round(spikes*1000);

sta=zeros(height,width,num_frames); %height, width, frames back
% stv=zeros(height,width,num_frames); %height, width, frames back
sta_store = zeros(height,width, num_frames, length(spikes), 3);

tic
icnt=0;

for i=spikes'
    
    start=find(frame_times>i,1)-num_frames;
    if(start>000)
        icnt=icnt+1
        for j=1:num_frames
            F = round(mvi.getFrame(start+j).getBuffer);
            sta(:,:,j) = sta(:,:,j) + round(reshape(F(1:3:end),width,height)'-0.5)+round(reshape(F(2:3:end),width,height)'-0.5)+round(reshape(F(3:3:end),width,height)'-0.5);
            sta_store(:,:,j, icnt,1) = double(round(reshape(F(1:3:end),width,height)'-0.5));
            sta_store(:,:,j, icnt,2)=double(round(reshape(F(2:3:end),width,height)'-0.5));
            sta_store(:,:,j, icnt,3) =double(round(reshape(F(3:3:end),width,height)'-0.5));
        end
    end
end
sta=sta/icnt;

% choose first frame to show
[junk,start_index] = max(sum(reshape(sta.^2,[],size(sta,3)),1));

% normalize STA color
sta = norm_image(sta);

% create slider control
%ha = make_loop_slider_list(start_index, 1, size(sta, 3), {@slider_plot, sta});

% plot once before any clicks
%slider_plot(ha, [], sta);
sta = repmat(sta, 1,1,1,3);
sta = permute(sta,[1 2 4 3]);
[sig_stixels] = significant_stixels(sta);
[row, col] = find(sig_stixels);
sig_stixels = [row,col];
% ranges = range(sta,3);
% [val, idx] = max(ranges(:));
% [x,y] = ind2sub(size(ranges), idx)
% x_iter = [x-1, x-1, x-1, x, x, x, x+1, x+1, x+1];
% y_iter = [y-1, y, y+1, y-1, y, y+1, y-1, y, y+1];
red = [];
blue = [];
green= [];
for iter = 1: size(sig_stixels,1)
    x_now = sig_stixels(iter,1);
    y_now = sig_stixels(iter,2);
    time = 0:-refresh:-refresh*(num_frames-1);
    sta_valueR = zeros(num_frames,size(sta_store, 4));
    sta_valueG = zeros(num_frames,size(sta_store, 4));
    sta_valueB = zeros(num_frames,size(sta_store, 4));
    
    for i = 1:num_frames
        for icnt = 1:size(sta_store, 4)
            sta_valueR(i, icnt) = sta_store(x_now,y_now,i,icnt,1);
            sta_valueG(i, icnt) = sta_store(x_now,y_now,i,icnt,2);
            sta_valueB(i, icnt) = sta_store(x_now,y_now,i,icnt,3);
        end
        i
    end
    red_temp = mean(sta_valueR,2)';
    green_temp = mean(sta_valueG,2)';
    blue_temp = mean(sta_valueB,2)';
    red = [red; red_temp];
    green = [green; green_temp];
    blue = [blue; blue_temp];
end
red_avg = mean(red,1);
green_avg = mean(green,1);
blue_avg = mean(blue,1);
colors = [red_avg; green_avg; blue_avg];


colors_flipped = fliplr(colors);
figure
plot(time, colors_flipped(1,:), 'r')
hold on
plot(time, colors_flipped(2,:), 'g')
plot(time, colors_flipped(3,:), 'b')

% this now gives the timecourse for each color for a different large cell
% Now analysis the cell you actually care about

%% STA calculation
%clearvars datarun color_super_large
%datarun.names.rrs_neurons_path='/Volumes/Analysis/2008-12-12-1/data006-nwpca/data006/data006.neurons';
%mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
target_cell2 = 3950;

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);
datarun=load_data(datarun,opt);

triggers=datarun.triggers; %onsets of the stimulus presentation

[mov,height,width,duration,refresh] = get_movie_ath(mdf_file,...
    triggers, 1,2);


[mvi] = load_movie(mdf_file, triggers);

%frames per trigger

bt_triggers = triggers - [0;triggers(1:end-1)];
avg_bt_triggers = mean(bt_triggers(2:end));
frames_per_trigger = round(avg_bt_triggers*1000/refresh);

last_trigger_time = ceil(triggers(end));

frame_times = zeros(ceil(last_trigger_time/refresh*1000),1);

% frame_times(1:25:end) = triggers

temp = linspace(triggers(1),refresh/1000,frames_per_trigger);
for i = 1: length(triggers)
    temp = linspace(triggers(i),triggers(i)+ (frames_per_trigger-1)*refresh/1000,frames_per_trigger)';
    frame_times(i*frames_per_trigger-(frames_per_trigger-1):i*frames_per_trigger) = temp;
end
frame_times = frame_times*1000; % in ms
cellID=find(datarun.cell_ids==target_cell2);


spikes=datarun.spikes{cellID};

% spikes are in s. convert to ms
spikes=round(spikes*1000);


sta=zeros(height,width,num_frames); %height, width, frames back
stv = zeros(height,width, num_frames, 3, length(spikes));

tic
icnt=0;

for i=spikes'
    
    start=find(frame_times>i,1)-num_frames;
    if(start>000)
        icnt=icnt+1
        for j=1:num_frames
            F = round(mvi.getFrame(start+j).getBuffer);
            % Change sta and stvto be -1 or +1
            sta(:,:,j) = sta(:,:,j) + round(reshape(F(1:3:end),width,height)'-0.5)+round(reshape(F(2:3:end),width,height)'-0.5)+round(reshape(F(3:3:end),width,height)'-0.5);
            stv(:,:,j, 1, icnt) = double(round(reshape(F(1:3:end),width,height)'-0.5));
            stv(:,:,j, 2, icnt)=double(round(reshape(F(2:3:end),width,height)'-0.5));
            stv(:,:,j, 3, icnt) =double(round(reshape(F(3:3:end),width,height)'-0.5));
        end
    end
end

%% Get significant stixels for histogram
sta=sta/icnt;
% normalize STA color
sta = norm_image(sta);
sta_rep = repmat(sta, 1,1,1,3);
sta_rep = permute(sta_rep,[1 2 4 3]);
[sig_stixels] = significant_stixels(sta_rep);
[row, col] = find(sig_stixels);
sig_stixels = [row,col];



%% Get timecourse in right dimensions
color_large = zeros(height, width, size(colors,2),3);
for i = 1:length(red_avg)
    color_large(:,:,i, 1) = repmat(colors(1,i), height, width);
end

for i = 1:length(green_avg)
    color_large(:,:,i, 2) = repmat(colors(2,i), height, width);
end

for i = 1:length(blue_avg)
    color_large(:,:,i,3) = repmat(colors(3,i), height, width);
end
color_super_large = repmat(color_large, 1,1,1,1,icnt);

%% Take inner product of timecourse and stv
temp_color = color_super_large .* stv(:,:,:,:, 1:icnt);
inner_prod_stv = squeeze(sum(temp_color,3));

%% Temporally weighted STV and skewness
stv_weight = zeros(height,width, 3);
stv_skew = zeros(height,width, 3);
stv_man = zeros(height,width, 3);

for i = 1:3
    stv_temp = squeeze(inner_prod_stv(:,:,i,:));
    for x = 1:height
        for y = 1:width
            stv_weight(x,y, i) = var(squeeze(inner_prod_stv(x,y,i,:)));
            stv_man(x,y, i) = sum(squeeze(inner_prod_stv(x,y,i,:)).*2)/icnt;%var(squeeze(inner_prod_stv(x,y,i,:)));
            stv_skew(x,y, i) = skewness(squeeze(inner_prod_stv(x,y,i,:)));
        end
    end
end
%% Normal STV, not scaled
stv_unscl = zeros(height, width, num_frames, 3);
for c = 1:3
    for i = 1:num_frames
        stv_temp = squeeze(stv(:,:,i,c,:));
        for x = 1:height
            for y = 1:width
                stv_unscl(x,y, i, c) = var(squeeze(stv_temp(x,y,:)));
            end
        end
    end
end


%% STV of each frame temporally scaled
stv_sep = zeros(height, width, num_frames, 3);
for color = 1:3
    for i = 1:num_frames
        for x = 1:height
            for y= 1:width
                stv_temp = squeeze(stv(x,y,i,color,:));
                stv_sep(x,y,i, color) = var(stv_temp); % look at all 30 frames not collasped but still temporally filtered
            end
        end
    end
end




%% Temporally scaled STA
sta = repmat(sta, 1,1,1,3);
temp_sta = color_large .* sta;
inner_prod_sta = squeeze(sum(temp_sta(:,:,:,2),3));


%% normal sta
figure
for j=1:num_frames
    imagesc(sta(:,:,j));
    colormap gray
    caxis([min(sta(:)),max(sta(:))]);
    colorbar
    title(['Normal STA Frame: ' num2str(j)]);
    pause(10/120);
    
end

%% temporally weighted sta
figure
imagesc(inner_prod_sta);
title('Temporally Weighted STA')
colormap gray
caxis([min(inner_prod_sta(:)),max(inner_prod_sta(:))]);
colorbar
pause(1/120);

%% unscaled stv, specify channel

figure;
for i = 1:num_frames
    imagesc(stv_unscl(:,:,i, 2));
    colormap gray
    stv_one_color = squeeze(stv_unscl(:,:,:,2));
    caxis([min(stv_one_color(:)),max(stv_one_color(:))]);
    colorbar
    title(['Normal STV Frame: ', num2str(i)])
    pause(10/120)
end



%% temporally weighted stv all frames collasped

figure;
imagesc(stv_weight(:,:,2));
title('Temporally Weighted STV')
colormap gray
caxis([min(stv_weight(:)),max(stv_weight(:))]);
colorbar

%% temporally weighted stv all frames collasped, mean = 0

figure;
imagesc(stv_man(:,:,2));
title('Temporally Weighted STV Deviation from 0')
colormap gray
caxis([min(stv_man(:)),max(stv_man(:))]);
colorbar


%% temporally weighted stv skewness

figure;
imagesc(stv_skew(:,:,2))
title('Temporally Weighted Skewness')
colormap gray
caxis([min(stv_skew(:)),max(stv_skew(:))]);
colorbar


%% temporally weighted stv looking at each frame individually
figure;
for i = 1:num_frames
    imagesc(stv_sep(:,:,i,2));
    colormap gray
    caxis([min(stv_sep(:)),max(stv_sep(:))]);
    colorbar
    title(['Temporally Weighted STV Frame: ' num2str(i)])
    pause(10/120)
end

%% histogram of peak spikes/ noise spikes at peak frame
figure
x_min = min(sig_stixels(:,1));
x_max = max(sig_stixels(:,1));
y_min = min(sig_stixels(:,2));
y_max = max(sig_stixels(:,2));
loc = 1;
for x = x_min:x_max
    for y = y_min:y_max
        if ismember([x,y], sig_stixels, 'rows');
            
            weight = 'bold';
        else
            weight = 'normal';
        end
    
        subplot(x_max-x_min+1,y_max-y_min+1, loc)
        hist(squeeze(inner_prod_stv(x,y,2, :)));
        title(sprintf('pixel (%d,%d)', x,y), 'fontsize', 8, 'fontweight', weight)
        set(gca, 'fontsize', 8)
        loc = loc+1;
    end
end





