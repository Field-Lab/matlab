% compute STA from CBP

% load spikes
clear
load ('/Volumes/Lab/Users/bhaishahster/NSEM_dataset_2005-04-26-0_cell1945_data009_share.mat')
spk_tm = spk_tm/20000; % seconds
spikes = spk_tm*1000; % in ms



datarun.names.rrs_neurons_path='/Volumes/Analysis/2005-04-26-0/data009-dum/data009/data009.neurons';
mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-20-1-0.48-33333.xml';
num_frames = 30; % both have to be run with the name number of frames
target_cell = 31;
% target_cell2 = 3169;


opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);
datarun=load_data(datarun,opt);

triggers=datarun.triggers; %onsets of the stimulus presentation
% [mov,height,width,duration,refresh] = get_movie_ath(mdf_file,...
%     triggers, 1,2);


[mvi] = load_movie(mdf_file, triggers);
height = double(mvi.getHeight);
width = double(mvi.getWidth);
duration = double(mvi.size);
refresh = double(mvi.getRefreshTime);

% refresh = 120;
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

% spikes=datarun.spikes{cellID};
% spikes are in s. convert to ms
% spikes=round(spikes*1000);

sta=zeros(height,width,3, num_frames); %height, width, frames back
% stv=zeros(height,width,num_frames); %height, width, frames back
sta_store = zeros(height,width, 3, num_frames, length(spikes));

tic
icnt=0;
length_spikes = length(spikes)

for i=spikes'
    
    start=find(frame_times>i,1)-num_frames;
    if(start>000)
        icnt=icnt+1;
        if mod(icnt, 1000) == 0
            disp(icnt)
        end
        
        for j=1:num_frames
            F = round(mvi.getFrame(start+j).getBuffer);
%             sta(:,:,j) = sta(:,:,j) + round(reshape(F(1:3:end),width,height)'-0.5)+round(reshape(F(2:3:end),width,height)'-0.5)+round(reshape(F(3:3:end),width,height)'-0.5);
                        sta(:,:,1, j) = sta(:,:,1,j) + round(reshape(F(1:3:end),width,height)'-0.5); % store the three color channels
                        sta(:,:,2, j) = sta(:,:,2,j) + round(reshape(F(2:3:end),width,height)'-0.5);
                        sta(:,:,3, j) = sta(:,:,3,j) + round(reshape(F(3:3:end),width,height)'-0.5);
            
            
            sta_store(:,:,1, j, icnt) = double(round(reshape(F(1:3:end),width,height)'-0.5));
            sta_store(:,:,2, j, icnt)=double(round(reshape(F(2:3:end),width,height)'-0.5));
            sta_store(:,:,3, j, icnt) =double(round(reshape(F(3:3:end),width,height)'-0.5));
        end
    end
end
sta=sta/icnt;

% normalize STA color
sta = norm_image(sta);

% choose first frame to show
[junk,start_index] = max(sum(reshape(sta.^2,[],size(sta,3)),1));

sta_store = sta_store(:,:,:,:,1:icnt);


% create slider control
%ha = make_loop_slider_list(start_index, 1, size(sta, 3), {@slider_plot, sta});

% plot once before any clicks
%slider_plot(ha, [], sta);
% sta = repmat(sta, 1,1,1,3);
% sta = permute(sta,[1 2 4 3]);
[sig_stixels] = significant_stixels(sta, 'select', 'thresh', 'thresh', 4);
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
        for icnt = 1:size(sta_store, 5)
            sta_valueR(i, icnt) = sta_store(x_now,y_now,1, i,icnt);
            sta_valueG(i, icnt) = sta_store(x_now,y_now,2,i,icnt);
            sta_valueB(i, icnt) = sta_store(x_now,y_now,3, i,icnt);
        end
        
    end
    red_temp = mean(sta_valueR,2)';
    green_temp = mean(sta_valueG,2)';
    blue_temp = mean(sta_valueB,2)';
    red = [red; red_temp];
    green = [green; green_temp];
    blue = [blue; blue_temp];
    iter
end
red_avg = mean(red,1);
green_avg = mean(green,1);
blue_avg = mean(blue,1);
colors = [red_avg; green_avg; blue_avg];

time = -(num_frames-1)*refresh:refresh:0
% colors_flipped = fliplr(colors);
figure
plot(time,colors(1,:), 'r')
hold on
plot(time,colors(2,:), 'g')
plot(time,colors(3,:), 'b')

%% normal sta
figure
% for j=1:num_frames
j =  23;
    imagesc(squeeze(sta(:,:,:,j)));
%     colormap gray
%     caxis([min(sta(:)),max(sta(:))]);
    colorbar
    title(['Normal STA Frame: ' num2str(j)]);
    pause(10/120);
    
% end


