clear
%% Get timecourse of related cell 
datarun.names.rrs_neurons_path='/Volumes/Analysis/2007-01-23-5/data001-map-data010/data001-map-data010.neurons';
mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-16-4-0.48-22222.xml';
num_frames = 32;
target_cell = 322;


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
time = 0:-refresh:-refresh*(num_frames-1);
sta_valueR = zeros(num_frames,size(sta_store, 4));
sta_valueG = zeros(num_frames,size(sta_store, 4));
sta_valueB = zeros(num_frames,size(sta_store, 4));

for i = 1:num_frames
    for icnt = 1:size(sta_store, 4)
    sta_valueR(i, icnt) = sta_store(x,y,i,icnt,1);
    sta_valueG(i, icnt) = sta_store(x,y,i,icnt,2);
    sta_valueB(i, icnt) = sta_store(x,y,i,icnt,3);
    end
    i
end
figure
plot(time, fliplr(mean(sta_valueR,2)'), 'r')
hold on
plot(time, fliplr(mean(sta_valueG,2)'), 'g')
plot(time, fliplr(mean(sta_valueB,2)'), 'b')

%% STA calculation

datarun.names.rrs_neurons_path='/Volumes/Analysis/2007-01-23-5/data010/data010-colleen/data010-colleen.neurons';
num_frames = 30;
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);
datarun=load_data(datarun,opt);

triggers=datarun.triggers; %onsets of the stimulus presentation

mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-8-4-0.48-11111.xml';
% triggers = [triggers; triggers + 1800];
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
    
% frame_times(1:24) = frame_times(1) + refresh/1000;
% figure
% imagesc(mov(:,:,2))

cellID=find(datarun.cell_ids==876)


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
        stv_store(:,:,j, icnt,1) = double(reshape(F(1:3:end),width,height)');
        stv_store(:,:,j, icnt,2)=double(reshape(F(2:3:end),width,height)');
        stv_store(:,:,j, icnt,3) =double(reshape(F(3:3:end),width,height)'); 
        end
    end
end
sta=sta/icnt;
% stv_store = stv_store(:,:,:,1:icnt, :);

frame_var = zeros(height, width, num_frames, 3);
for color = 1:3
for frame = 1:num_frames
    
    stv_frame = squeeze(stv_store(:,:,frame,1:icnt, color));
    for x = 1:height
        for y = 1:width
            frame_var(x,y,frame, color) = var(squeeze(stv_frame(x,y,:)));
        end
    end
end
color
end



% See STA .. STA should be good and its a proof that the code is working
% fine ..
figure;
for j=26%1:num_frames

imagesc(sta(:,:,j));
colormap gray
caxis([min(sta(:)),max(sta(:))]);
colorbar
pause(10/120);
end

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
