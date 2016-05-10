
clear
tic
%% ------------------------------ INPUTS -----------------------------------
cells = {3,1147,1188,2208,2643,2866,3512,4863,5060,5915,6706,7636}; % from vision
file_name = '2015-11-09-8/data006/data006'; % classification file
screen_width = 640; % in pixels 
screen_height = 320; 
stixels_ref = 8; % stixel size of white noise run
stixels_mask = 1; % = 2 for NSEM
raster_length = 15;
LES = 1;
mask_type = 'sig'; % options are sig: only significant stixels, or est: approx from STA fit
filter_type = 'sig'; % same options as above
params.thresh = 3;

%% ------------------------------- Save details of movie in a new folder ------------------------------------------
folder = ['/Volumes/Lab/Users/Nora/new_stim_nora/mask_NSEM/mask_movie_' num2str(cells{1}) '_11_stix' num2str(stixels_mask) '/'];
mkdir(folder); 
save([folder '/inputs.mat'])

%% ------------------------------- Load Data ------------------------------------------
datarun.names.rrs_params_path=['/Volumes/Acquisition/Analysis/', file_name, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Acquisition/Analysis/', file_name, '.sta'];
opt=struct('verbose',0,'load_params',1,'load_neurons',0,'load_obvius_sta_fits',false, 'load_sta', 1, 'load_sta_params', 0, 'load_all',false);
opt.load_sta_params.frames = 1:30;% if this line is missing, will error; have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);

%% ------------------------------- Plot Vision STA -----------------------------
mask_width = screen_width/stixels_mask; 
mask_height = screen_height/stixels_mask;
mask = zeros(mask_height, mask_width, length(cells));
linear_filter = zeros(mask_height, mask_width, length(cells));

for i_cell = 1:length(cells)
    myMap = zeros(mask_height, mask_width); % pixesl on the screen
    [cell_numbers, cell_type, cell_type_number] = get_cell_indices(datarun, cells{i_cell});
    
    the_fit = datarun.stas.fits{cell_numbers};
    ctr = the_fit.mean;
    rad = mean(the_fit.sd);
    % axis([0 screen_width/stixels_ref 0 screen_height/stixels_ref])
    if strcmp(mask_type, 'est')
        [X,Y] = drawEllipse_upsampled([ctr [rad rad] the_fit.angle]);
        X_large =  round(X*stixels_ref/stixels_mask);
        Y_large =  round(Y*stixels_ref/stixels_mask);
        for i = 1:length(X_large)
            myMap(Y_large(i),X_large(i)) = 1;
        end
        mask(:,:,i_cell) = imfill(myMap,'holes');
    elseif strcmp(mask_type, 'sig')
        A = significant_stixels(datarun.stas.stas{cell_numbers}, params);
        A = ExtractNLargestBlobs(full(A),1);
        mask(:,:,i_cell) = imresize(full(A), stixels_ref/stixels_mask, 'nearest');
    end

    % Make the linear filter
    if LES
        mkdir([folder 'LES/']);
        if strcmp(filter_type, 'est')
            ctr = round(stixels_ref/stixels_mask*the_fit.mean);
            rad = stixels_ref/stixels_mask*mean(the_fit.sd);
            filter = zeros(mask_height, mask_width);
            filter_size = ceil(3*rad)+1-mod(ceil(3*rad),2); % weird other stuff makes it odd
            alpha = (filter_size-1)/(2*rad); % just from matlab's weird code definition for gausswin
            idx = -floor(filter_size/2):floor(filter_size/2);
            filter(idx + ctr(2), idx + ctr(1)) = gausswin(filter_size, alpha)*gausswin(filter_size, alpha)';
            filter(~mask(:,:,i_cell)) = 0;
            linear_filter(:,:,i_cell) = filter;
        elseif strcmp(filter_type, 'sig')
            [A,~,filter] = significant_stixels(datarun.stas.stas{cell_numbers}, params);
            A = ExtractNLargestBlobs(full(A),1);
            filter(~A) = 0; % outside significant stix = 0;
            filter = imresize(filter, stixels_ref/stixels_mask, 'nearest');
            filter(~logical(mask(:,:,i_cell))) = 0; % outside mask is also 0
            linear_filter(:,:,i_cell) = filter;
        end
        % figure; imagesc(linear_filter(:,:,i_cell)); axis image
        %L1 = sum(linear_filter(:));
    end
    
end

% cut movie down for single stixel
if stixels_mask == 1
    idx_x = (1:mask_height/2) + mask_height/4;
    idx_y = (1:mask_width/2) + mask_width/4;
    check_mask = zeros(mask_height, mask_width);
    check_mask(idx_x, idx_y) = 1;
    total_mask = sum(mask, 3);
    if any(total_mask(~check_mask))
        warning('Some RF is being cut off')
    end
    mask = mask(idx_x, idx_y, :);
    linear_filter = linear_filter(idx_x, idx_y, :);
end
figure; imagesc(sum(mask,3))



%% ------------------------------- Create original NSEM movie -----------------------------

% Load up full movie
i = 1;
total_frames = raster_length*120
current_frame = 0;
while current_frame < total_frames;
    load(['/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSinterval/matfiles/movie_chunk_' num2str(i) '.mat'])
    if current_frame + size(movie_chunk,3) > total_frames
        movie_chunk = movie_chunk(:,:,1:(total_frames - current_frame));
    end
    current_frame = current_frame + size(movie_chunk,3);
    movie_chunk = movie_chunk.*repmat(sum(mask,3)', [1,1,size(movie_chunk,3)]);
    movie_chunk(~repmat(sum(mask,3)', [1,1,size(movie_chunk,3)])) = 64;
    save([folder 'movie_chunk_' num2str(i) '.mat'], 'movie_chunk')
    if LES
        for i_frame = 1:size(movie_chunk,3)
            frame = movie_chunk(:,:,i_frame);
            for i_cell = 1:length(cells)
                drive = sum(sum(linear_filter(:,:,i_cell)' .* frame));
                LES_value = round(drive/sum(sum(linear_filter(:,:,i_cell))));
                frame(logical(mask(:,:,i_cell)')) = LES_value;
            end
            movie_chunk(:,:,i_frame) = frame;
        end
        save([folder 'LES/movie_chunk_' num2str(i) '.mat'], 'movie_chunk')
    end
    i = i+1;
end

moviename = ['movie_' num2str(cells{1}) '.rawMovie'];
generate_maskmovie_from_mat(folder, moviename, current_frame)
%orig_stim = get_rawmovie([folder moviename], 240, 0);
if LES
    moviename = ['movie_' num2str(cells{1}) '_LES.rawMovie'];
    generate_maskmovie_from_mat([folder 'LES/'], moviename, current_frame)
    %LES_stim = get_rawmovie([[folder 'LES/'] moviename], 240, 0);
end

%{
for i = 1:size(orig_stim,1)
    subplot(2,1,1)
    imagesc(squeeze(orig_stim(i,:,:))')
    caxis([0 255])
    colormap gray
    axis image
    if LES
        subplot(2,1,2)
        imagesc(squeeze(LES_stim(i,:,:))')
        caxis([0 255])
        colormap gray
        axis image
    end
    pause(0.008)
    
end

drive_1 = sum(sum(squeeze(double(orig_stim(end,:,:))).*linear_filter(:,:,1)'))
drive_2 = sum(sum(squeeze(double(LES_stim(end,:,:))).*linear_filter(:,:,1)'))

drive_3 = sum(sum(squeeze(double(orig_stim(end,:,:))).*linear_filter(:,:,2)'))
drive_4 = sum(sum(squeeze(double(LES_stim(end,:,:))).*linear_filter(:,:,2)'))
%}
toc

%% ------------------------------- Calculate LES -----------------------------






