clear 

dataparam.date='2016-02-17-6/';
dataparam.concatname='data014';
% dataparam.mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-1-1-0.16-11111-11x1-119.5.xml';
% dataparam.mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-gaussian-1-1-0.16-11111-14x1-119.5.xml';
% dataparam.mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-gaussian-1-1-0.16-11111-97x1-119.5.xml';
dataparam.mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-gaussian-1-1-0.16-11111-14x1-119.5.xml';

vision_id = 1928;

image_width = 14;
image_height = 1;
num_frames = 30;
num_colors=3;
dataparam.file_name_right = [dataparam.date, '/', dataparam.concatname,'/', dataparam.concatname];


datarun.names.rrs_neurons_path=['/Volumes/Analysis/', dataparam.file_name_right, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', dataparam.file_name_right, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Analysis/', dataparam.file_name_right, '.sta'];

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false, 'load_sta', 0, 'load_sta_params', 1, 'load_all',true);
opt.load_sta_params.save_rf = 1;
% opt.load_sta_params.frames = 1:fitparam.num_frames;% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);

triggers = datarun.triggers;
cell_indices = get_cell_indices(datarun, vision_id);
spikes = {datarun.spikes{cell_indices}};

% datarun.triggers= datarun.triggers(1:2);
[inputs, refresh, duration] = get_wn_movie_ath(datarun, dataparam.mdf_file, 0);

real_frame = zeros(image_width, image_height, num_colors);
real_frame(:,:,1,1) = reshape(inputs(1:3:image_width*image_height*3)',image_width, image_height);
real_frame(:,:,2,1) = reshape(inputs(2:3:image_width*image_height*3)',image_width, image_height);
real_frame(:,:,3,1) = reshape(inputs(3:3:image_width*image_height*3)',image_width, image_height);

disp('done processing movie')
pointer = image_width*image_height*num_colors+1;
%     pointer = 2+25+2;
i =2;
while pointer+image_height*image_width*3-1<size(inputs,2)*size(inputs,1)
    temp = inputs(pointer:pointer+image_height*image_width*3-1);
    real_frame(:,:,1,i) = reshape(temp(1:3:end), image_width, image_height);
    real_frame(:,:,2,i) = reshape(temp(2:3:end), image_width, image_height);
    real_frame(:,:,3,i) = reshape(temp(3:3:end), image_width, image_height);
    
    pointer = pointer+image_height*image_width*3;
    i = i+1;
end


length_of_time = ceil(triggers(end))+1;
upsampled_num_frames = length_of_time*120;

upsample_factor = round(refresh/(100/12));
bt_triggers = triggers - [0;triggers(1:end-1)];
avg_bt_triggers = mean(bt_triggers(2:end));
triggers = [triggers; triggers(end) + avg_bt_triggers];

for j = 1:size(spikes,2)
    
    
    frames_needed{j} = zeros(3,(length(triggers)-2)*100+120);
    a = kron(1:upsampled_num_frames/upsample_factor, ones(1,upsample_factor));
    frames_needed{j}(1,:) = a(1:(length(triggers)-2)*100+120);
    
    %     frames_needed{j}(1,:) = kron(1:upsampled_num_frames/upsample_factor, ones(1,upsample_factor));
    
    
    for i= 1:length(triggers)-1
        spacing = linspace(triggers(i), triggers(i+1),101);
        frames_needed{j}(2, (i-1)*100+1:(i-1)*100+100)= spacing(1:end-1);
    end
    
    for i = 1:size(frames_needed{j},2)-1
        frames_needed{j}(3,i) = sum(spikes{j} >= frames_needed{j}(2,i) & spikes{j} < frames_needed{j}(2,i+1));
    end
end


sta = cell(size(spikes,2),1);
for i = 1:size(spikes,2)
    sta{i} =zeros(size(real_frame,1),size(real_frame,2),num_colors, num_frames);
end



for cel = 1:size(spikes,2)
    
    for i = 1:size(frames_needed{cel},2)
        %         if start_points(j) -1+ i <= length(spikes_by_frame)
        if frames_needed{cel}(3, i) ~= 0
            
            if i > num_frames
                
                for t = 1:num_frames
                    subtract = num_frames - t +1;
                    sta{cel}(:,:, :,subtract) = sta{cel}(:,:, :, subtract) + real_frame(:,:,:,i-t) * frames_needed{cel}(3, i);
                end
            end
        end
    end
    
end


%
% % fitmovie = real_frame;
% tstim = refresh/120;
%
% movie_size = size(real_frame);
% STA = zeros(movie_size(1),movie_size(2),30);
% fitframes = movie_size(4);
% % STA = [];
%
% for i = 1:length(fitspikes)
%     sp_frame = floor(fitspikes(i)/tstim);
%     if sp_frame > 29 && sp_frame<fitframes
%         STA = STA+double(real_frame(:,:,(sp_frame-29):sp_frame));
%     end
% end
% STA = STA./length(fitspikes);
% for i = 1:30
%    imagesc(STA(:,:,i)')
%    colormap gray
%    axis image
%    title('You should see an STA here')
%    pause(0.1)
% end
%
% STA_full = STA;
%
% STA = squeeze(STA);
% STA = abs(sum(STA(:,:,24:27), 3));
% row = max(STA);
% col = max(STA');
% x = find(row == max(row));
% y = find(col == max(col));
% center = [y x];
%
% imagesc(STA')
% axis image
%
% if center_verification
%     title('Click on the center of the STA')
%     [x, y] = ginput(1);
%     center = round([x y]);
% else
%     hold on
%     plot(center(1), center(2), 'r*')
%     title('The red dot should be over the center of the STA')
%     pause(0.1)
% end
%
%
% end

figure; plot(squeeze(sta{1}(3,1,2,:)))

