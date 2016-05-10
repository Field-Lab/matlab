clear
tic
%% ------------------------------ INPUTS -----------------------------------
cells = {3,1147,1188,2208,2643,2866,3512,4863,5060,5915,6706,7636}; % from vision
piece = '2015-11-09-8';
run = '006'; % classification file
screen_width = 640; % in pixels 
screen_height = 320; 
stixels_ref = 8; % stixel size of white noise run
stixels_mask = 8; % = 2 for NSEM
raster_length = 15;
LES = 1;
mask_type = 'sig'; % options are sig: only significant stixels, or est: approx from STA fit
filter_type = 'sig'; % same options as above
params.thresh = 3;

%% ------------------------------- Load Data ------------------------------------------
datarun =load_data([piece '/streamed/data' run '/data' run]);
datarun = load_sta(datarun);
datarun = load_params(datarun);


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

%%
mask = sum(mask,3);
figure; imagesc(mask)

%% average STAs
off_avg = compute_average_rf(datarun, 'Off Parasol');
on_avg = compute_average_rf(datarun, 'On Parasol');
figure; imagesc(sum(off_avg,3)); title('OFF');
figure; imagesc(sum(on_avg,3)); title('ON')

%%
for i_cell = 1:length(cells)
    [cell_numbers, cell_type, cell_type_number] = get_cell_indices(datarun, cells{i_cell});
    the_fit = datarun.stas.fits{cell_numbers};
    ctr = round(stixels_ref/stixels_mask*the_fit.mean);
    idx1 = (ctr(2)-7):(ctr(2)+7);
    idx2 = (ctr(1)-7):(ctr(1)+7);
    centered_mask = mask();
    figure(7); imagesc(imresize(centered_mask, 8));
    pause()
end


