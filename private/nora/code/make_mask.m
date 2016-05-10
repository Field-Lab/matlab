function final_mask = make_mask(cells, classification_file, varargin)


p = inputParser;
p.addParameter('screen_size', [640, 320])
p.addParameter('stixel_class', 8)
p.addParameter('stixel_mask', 2)
p.addParameter('mask_type', 'est') % options are sig: only significant stixels, or est: approx from STA fit
p.addParameter('sigmas',  2);
p.addParameter('circular', 0)

% The rest of these only matter for the linearly equivalent stimulus, so
% you can ignore if you don't want that
p.addParameter('LES_movie', 'NSinterval') % the movie to compare
p.addParameter('LES_length', 0) % the length of the LES
p.addParameter('filter_type', 'est') % options are sig: only significant stixels, or est: approx from STA fit
p.addParameter('thresh', 3)
p.addParameter('save_folder', '/Volumes/Lab/Users/Nora/new_stim_nora/mask_NSEM/test')
p.addParameter('complement', 1) % make complement with LES?
p.parse(varargin{:});

screen_width = p.Results.screen_size(1);
screen_height = p.Results.screen_size(2);
stixels_ref = p.Results.stixel_class; % stixel size of white noise run
stixels_mask = p.Results.stixel_mask; % = 2 for NSEM
raster_length = p.Results.LES_length;
if raster_length>0; LES = 1; else LES = 0; end
mask_type = p.Results.mask_type; % options are sig: only significant stixels, or est: approx from STA fit
filter_type = p.Results.filter_type; % same options as above
params.thresh = p.Results.thresh;
save_folder = p.Results.save_folder;

if ~strcmp(p.Results.LES_movie, {'NSbrownian', 'NSinterval'});
    warning('Linear equivalent stimulus is only set up for certain movies')
    LES = 0;
end


%% ------------------------------- Save details of movie in a new folder ------------------------------------------
if LES
    folder = [save_folder 'mask_' num2str(cells{1}) '_sigma' num2str(p.Results.sigmas) '/'];
    mkdir(folder);
    save([folder '/inputs.mat'])
    mkdir([folder 'LES/']);
    if p.Results.complement;  mkdir([folder 'comp_LES/']); end
end

%% ------------------------------- Load Data ------------------------------------------
datarun.names.rrs_params_path=['/Volumes/Acquisition/Analysis/', classification_file, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Acquisition/Analysis/', classification_file, '.sta'];
opt=struct('verbose',0,'load_params',1,'load_neurons',0,'load_obvius_sta_fits',false, 'load_sta', 1, 'load_sta_params', 0, 'load_all',false);
opt.load_sta_params.frames = 1:6;% if this line is missing, will error; have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);

%% ------------------------------- Plot Vision STA -----------------------------
mask_width = screen_width/stixels_mask;
mask_height = screen_height/stixels_mask;
mask = zeros(mask_height, mask_width, length(cells));
linear_filter = zeros(mask_height, mask_width, length(cells));

for i_cell = 1:length(cells)
    myMap = zeros(mask_height, mask_width); % pixesl on the screen
    [cell_numbers, ~, ~] = get_cell_indices(datarun, cells{i_cell});
    
    % axis([0 screen_width/stixels_ref 0 screen_height/stixels_ref])
    if strcmp(p.Results.mask_type, 'est')
        the_fit = datarun.stas.fits{cell_numbers};
        ctr = the_fit.mean;
        if p.Results.circular
            rad = p.Results.sigmas*geomean(the_fit.sd);
            [X,Y] = drawEllipse_upsampled([ctr [rad rad] the_fit.angle]);
        else
            rad = p.Results.sigmas*the_fit.sd;
            [X,Y] = drawEllipse_upsampled([ctr rad the_fit.angle]);
        end
        X_large =  round(X*stixels_ref/stixels_mask);
        Y_large =  round(Y*stixels_ref/stixels_mask);
        X_large(X_large > mask_width) = mask_width;
        Y_large(Y_large > mask_height) = mask_height;
        X_large(X_large < 1) = 1;
        Y_large(Y_large < 1) = 1;
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
        if strcmp(filter_type, 'est')
            ctr = round(stixels_ref/stixels_mask*the_fit.mean);
            rad = stixels_ref/stixels_mask*mean(the_fit.sd);
            filter = zeros(mask_height, mask_width);
            filter_size = ceil(5*rad)+1-mod(ceil(5*rad),2); % weird other stuff makes it odd
            alpha = (filter_size-1)/(2*rad); % just from matlab's weird code definition for gausswin
            idx = -floor(filter_size/2):floor(filter_size/2);
            idx_x = idx + ctr(2);
            idx_y = idx + ctr(1);
            idx_x = idx_x(idx_x >0);
            idx_y = idx_y(idx_y >0);
            idx_x = idx_x(idx_x <= mask_height);
            idx_y = idx_y(idx_y <= mask_width);
            filter(idx_x, idx_y) = gausswin(filter_size, alpha)*gausswin(filter_size, alpha)';
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
% if stixels_mask == 1
%     idx_x = (1:mask_height/2) + mask_height/4;
%     idx_y = (1:mask_width/2) + mask_width/4;
%     check_mask = zeros(mask_height, mask_width);
%     check_mask(idx_x, idx_y) = 1;
%     total_mask = sum(mask, 3);
%     if any(total_mask(~check_mask))
%         warning('Some RF is being cut off')
%     end
%     mask = mask(idx_x, idx_y, :);
%     linear_filter = linear_filter(idx_x, idx_y, :);
% end
final_mask = sum(mask,3);
final_mask(final_mask > 1) = 1;
figure; imagesc(final_mask); axis image



%% ------------------------------- Create LES movie -----------------------------

if LES
    % Load up full movie
    i = 1;
    total_frames = raster_length*120;
    current_frame = 0;
    while current_frame < total_frames;
        
        % Load up original movie
        load(['/Volumes/Data/Stimuli/movies/eye-movement/current_movies/' p.Results.LES_movie '/matfiles/movie_chunk_' num2str(i) '.mat'])
        if current_frame + size(movie_chunk,3) > total_frames
            movie_chunk = movie_chunk(:,:,1:(total_frames - current_frame));
        end
        current_frame = current_frame + size(movie_chunk,3);
        
        % Calculate linear drive for each frame for each cell
        for i_frame = 1:size(movie_chunk,3)
            frame = movie_chunk(:,:,i_frame);
            for i_cell = 1:length(cells)
                drive = sum(sum(linear_filter(:,:,i_cell)' .* frame));
                LES_value = round(drive/sum(sum(linear_filter(:,:,i_cell))));
                frame(logical(mask(:,:,i_cell)')) = LES_value;
            end
            movie_chunk(:,:,i_frame) = frame;
        end
        
        % save new matfiles
        if p.Results.complement; save([folder 'comp_LES/movie_chunk_' num2str(i) '.mat'], 'movie_chunk'); end
        movie_chunk(repmat(~logical(final_mask'), [1 1 size(movie_chunk,3)])) = 64;
        save([folder 'LES/movie_chunk_' num2str(i) '.mat'], 'movie_chunk')
        i = i+1;
    end
    
    % create LES rawmovie
    moviename = ['movie_' num2str(cells{1}) '_LES.rawMovie'];
    generate_maskmovie_from_mat([folder 'LES/'], moviename, current_frame)
    if p.Results.complement
        moviename = ['movie_' num2str(cells{1}) '_comp_LES.rawMovie'];
        generate_maskmovie_from_mat([folder 'comp_LES/'], moviename, current_frame)
    end
    save([folder '/mask_and_filter.mat'],'final_mask', 'mask', 'linear_filter')
end

end




