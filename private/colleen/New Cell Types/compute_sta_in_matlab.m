clear
%% Get timecourse of related cell 
datarun.names.rrs_neurons_path='/Volumes/Analysis/2007-01-23-5/data001-map-data010/data001-map-data010.neurons';
mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-16-4-0.48-22222.xml';
num_frames = 32;
target_cell = 8;
target_cell2 = 19;


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
        sta(:,:,j) = sta(:,:,j) + reshape(F(1:3:end),width,height)'+reshape(F(2:3:end),width,height)'+reshape(F(3:3:end),width,height)';
        sta_store(:,:,j, icnt,1) = double(reshape(F(1:3:end),width,height)');
        sta_store(:,:,j, icnt,2)=double(reshape(F(2:3:end),width,height)');
        sta_store(:,:,j, icnt,3) =double(reshape(F(3:3:end),width,height)'); 
        end
    end
end
sta=sta/icnt;

% choose first frame to show
[junk,start_index] = max(sum(reshape(sta.^2,[],size(sta,3)),1));

% normalize STA color
sta = norm_image(sta);

% create slider control
ha = make_loop_slider_list(start_index, 1, size(sta, 3), {@slider_plot, sta});

% plot once before any clicks
slider_plot(ha, [], sta);

ranges = range(sta,3);
[val, idx] = max(ranges(:));
[x,y] = ind2sub(size(ranges), idx)
x_iter = [x-1, x-1, x-1, x, x, x, x+1, x+1, x+1];
y_iter = [y-1, y, y+1, y-1, y, y+1, y-1, y, y+1];
red = [];
blue = [];
green= [];
for iter = 1: length(x_iter)
    x_now = x_iter(iter);
    y_now = y_iter(iter);
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
red_avg = (red_avg(:) - min(red_avg(:))) / ( max(red_avg(:)) - min(red_avg(:)));
green_avg = (green_avg(:) - min(green_avg(:))) / ( max(green_avg(:)) - min(green_avg(:)));
blue_avg = (blue_avg(:) - min(blue_avg(:))) / ( max(blue_avg(:)) - min(blue_avg(:)));

figure
plot(time, fliplr(red_avg), 'r')
hold on
plot(time, fliplr(green_avg), 'g')
plot(time, fliplr(blue_avg), 'b')

% this now gives the timecourse for each color for a different large cell
% Now analysis the cell you actually care about

%% STA calculation

cellID=find(datarun.cell_ids==target_cell2);


spikes=datarun.spikes{cellID};
% spikes = spikes(1:20000);
% spikes are in s. convert to ms
spikes=round(spikes*1000);
%fr - frames - are in ms
%fr=round(triggers(1)*1000:refresh:triggers(end)*1000);
% make fr better
% refresh = 33.309894493290145;
% fr=[];
% for itrig=1:length(triggers)
% fr=[fr,triggers(itrig)*1000+refresh];
% end
% fr=fr';

sta=zeros(height,width,num_frames); %height, width, frames back
stv=zeros(height,width,num_frames); %height, width, frames back
stv_store = zeros(height,width, num_frames, length(spikes), 3);

tic
icnt=0;

for i=spikes'
 
    start=find(frame_times>i,1)-num_frames; 
    if(start>000)
    icnt=icnt+1
        for j=1:num_frames
        F = round(mvi.getFrame(start+j).getBuffer);
        sta(:,:,j) = sta(:,:,j) + reshape(F(1:3:end),width,height)'+reshape(F(2:3:end),width,height)'+reshape(F(3:3:end),width,height)';
        stv_store(:,:,j, icnt,1) = double(round(reshape(F(1:3:end),width,height)'-0.5));
        stv_store(:,:,j, icnt,2)=double(round(reshape(F(2:3:end),width,height)'-0.5));
        stv_store(:,:,j, icnt,3) =double(round(reshape(F(3:3:end),width,height)'-0.5)); 
        end
    end
end

red_large = zeros(height, width, length(red_avg));
for i = 1:length(red_avg)
    red_large(:,:,i) = repmat(red_avg(i), height, width);
end
red_super_large = repmat(red_large, 1,1,1,icnt);

temp_red = red_super_large .* squeeze(stv_store(:,:,:,1:icnt,1));

green_large = zeros(height, width, length(green_avg));
for i = 1:length(green_avg)
    green_large(:,:,i) = repmat(green_avg(i), height, width);
end
green_super_large = repmat(green_large, 1,1,1,icnt);

temp_green = green_super_large .* squeeze(stv_store(:,:,:,1:icnt,2));

blue_large = zeros(height, width, length(blue_avg));
for i = 1:length(blue_avg)
    blue_large(:,:,i) = repmat(blue_avg(i), height, width);
end
blue_super_large = repmat(blue_large, 1,1,1,icnt);

temp_blue = blue_super_large .* squeeze(stv_store(:,:,:,1:icnt,3));

inner_prod_red = squeeze(sum(temp_red,3));
inner_prod_green = squeeze(sum(temp_green,3));
inner_prod_blue = squeeze(sum(temp_blue,3));
inner_prod = zeros(height, width, icnt, 3);
inner_prod(:,:,:,1) = inner_prod_red;
inner_prod(:,:,:,2) = inner_prod_green;
inner_prod(:,:,:,3) = inner_prod_blue;

% stv_store = stv_store(:,:,:,1:icnt, :);

for color = 1:3

    for x = 1:height
        for y = 1:width
            frame_var(x,y, color) = var(squeeze(inner_prod(x,y,:, color)));
        end
    end
color
end


% green = temp_green(1,1,:,:);
% green = squeeze(green);
% mean_green = mean(green,2);
% figure; plot(mean_green);
% See STA .. STA should be good and its a proof that the code is working
% fine ..
% c1423 = squeeze(green_large(14,23,:));
% nc1423 = (c1423(:) - min(c1423(:))) / ( max(c1423(:)) - min(c1423(:)));
% c11 = squeeze(green_large(1,1,:));
% nc11 = (c11(:) - min(c11(:))) / ( max(c11(:)) - min(c11(:)));

% green = temp_green(14,23,:,:);
% green = squeeze(green);
% sum1423=repmat(nc1423,1, 23900).*green;
% green = temp_green(1,1,:,:);
% green = squeeze(green);
% sum11 = repmat(nc11,1, 23900).*green;
% figure;
% sum_green = var(green,0,2);

% temporally weighted variance 
figure; imagesc(frame_var(:,:,2)); colorbar


for j=26%1:num_framesfor j=26%1:num_frames

imagesc(sta(:,:,j));
colormap gray
caxis([min(sta(:)),max(sta(:))]);
colorbar
pause(1/120);
end

imagesc(sta(:,:,j));
colormap gray
caxis([min(sta(:)),max(sta(:))]);
colorbar
pause(1/120);


figure;
for j=2%1:num_frames
% rgbImage = cat(3, squeeze(frame_var(:,:,j,1)), squeeze(frame_var(:,:,j,2)), squeeze(frame_var(:,:,j,3)));
% imagesc(rgb2gray(rgbImage));
imagesc(frame_var(:,:,j,2))
colormap gray
caxis([min(frame_var(:)),max(frame_var(:))]);
colorbar
pause(10/120);
end


figure
time = 0:-refresh:-refresh*(num_frames-1);
ranges = range(sta,3);
[val, idx] = max(ranges(:));
[x,y] = ind2sub(size(ranges), idx)
% x= 12
% y = 16
x = 21;
y = 29;

sta_value = zeros(num_frames,1);
for i = 1:30
    sta_value(i) = sta(x,y,i);
end
figure
plot(time, fliplr(sta_value'))


time = 0:-refresh:-refresh*(num_frames-1);
x = 7;
y = 21;
stv_value = zeros(num_frames,1);
for i = 1:num_frames
    stv_value(i) = frame_var(x,y,i,2);
end
figure
plot(time, fliplr(stv_value'))

stimuli = zeros(icnt, 3);
bins = [0 0 0 ; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1];

for i = 1:icnt
    stimuli(i,:) = stv_store(x,y, 19, i, :);
end

for b = 1:size(bins,1)
    which_bin{b} = find(ismember(stimuli, bins(b,:), 'rows'))
end

% to_hist = 1:size(bins,1)';
for i = 1:size(bins,1)
    to_hist(1,i) = size(which_bin{i},1);
end

figure
for i =1:size(bins,1)
plot(i, to_hist(i)/icnt,  'ko-','markerfacecolor', bins(i,:))
hold on
end

set(gca, 'xticklabel' ,{'000', '001', '010', '011', '100', '101', '110', '111'}) 
title({'Max STA frame stimulus characteristics'; 'Cell 3946 Data 2007-03-27-2/data009/data009-colleen'})
ylabel('Proportion of frames of a particular color')
xlabel('RGB Color')
