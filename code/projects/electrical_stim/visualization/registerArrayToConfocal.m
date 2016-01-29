%% This script registers two images using control point selection
% LG 10/2015

experiment_id = '2015-10-06-6'; 

switch experiment_id
    case '2016-01-05-4'
        alignment_image = imread('/Volumes/Transfer/2015-01-05-4/Imaging/vasculature_stitched.tif'); 
        % Choose the channel containing vasculature only.
        vasculature_array = alignment_image(:,:,1);
        % Load the confocal vasculature image
        vasculature_confocal = imread('/Volumes/Transfer/2015-01-05-4/Imaging/confocal/vasc-2016-01-05-4.tif');
    case '2015-10-06-6'
        % Load the alignment image containing the array and the vasculature
        alignment_image = imread('/Volumes/Data/2015-10-06-6/Imaging/vasculature_alignment_with_array/2015-10-06-6_stitched_sm.jpg (blue).jpg');

        % Choose the channel containing vasculature only.
        vasculature_array = alignment_image(:,:,3);
        % Load the confocal vasculature image
        vasculature_confocal = imread('/Volumes/Data/2015-10-06-6/Imaging/confocal_vasculature.tif');
        confocal = imread('/Volumes/Data/2015-10-06-6/Imaging/confocal/2015-10-6-6_tiled_MIP_RGB.tiff');
        vasculature_confocal = vasculature_confocal(1:size(confocal,1),1:size(confocal,2));
        load('/Volumes/Analysis/2015-10-06-6/image analysis/coregistration_control_pts.mat');
        movingPoints = movingPoints_2015_10_06_6;
        fixedPoints = fixedPoints_2015_10_06_6;
        load('/Volumes/Analysis/2015-10-06-6/image analysis/electrode_XY_coords.mat'); 
        pna = confocal(:,:,1); 
        tubulin = confocal(:,:,2); 
        dapi = confocal(:,:,3); 
    case '2015-10-06-3'
        % Load the alignment image containing the array and the vasculature
        alignment_image = imread('/Volumes/Data/2015-10-06-3/Imaging/vasculature_alignment_with_array/2015-10-06-3_vasculature_stitch_light.tif');
       
        vasculature_array = rot90(alignment_image(:,:,1),2); clear alignment_image; 
        vasculature_confocal = imread('/Volumes/Data/2015-10-06-3/Imaging/confocal/MIP_tiled_hyperstack-pna.tif');
        tubulin = imread('/Volumes/Data/2015-10-06-3/Imaging/confocal/MIP_tiled_hyperstack-tubulin.tif');
        dapi = imread('/Volumes/Data/2015-10-06-3/Imaging/confocal/MIP_tiled_hyperstack-dapi.tif'); 
        pna = vasculature_confocal ; 
        load('/Volumes/Analysis/2015-10-06-3/image analysis/coregistration_control_points.mat');
    case '2015-11-09-3'
        alignment_image =imread('/Volumes/Data/2015-11-09-3/Imaging/2015-11-09-3-brightened.jpg');
        vasculature_array = alignment_image(:,:,1); clear alignment_image;
        vasculature_confocal = imread('/Volumes/Data/2015-11-09-3/Imaging/confocal/C3-2015-11-09-3_tilestack_MIP.tif');
        tubulin = imread('/Volumes/Data/2015-11-09-3/Imaging/confocal/C1-2015-11-09-3_tilestack_MIP.tif');
        dapi = imread('/Volumes/Data/2015-11-09-3/Imaging/confocal/C2-2015-11-09-3_tilestack_MIP.tif');
        pna = vasculature_confocal;
        load('/Volumes/Analysis/2015-11-09-3/image analysis/coregistration_control_pts.mat');
        load('/Volumes/Analysis/2015-11-09-3/image analysis/registered_electrode_coordinates.mat'); 
    case '2015-11-09-10'
        alignment_image = imread('/Volumes/Data/2015-11-09-10/Imaging/2015-11-09-10.jpg');
        vasculature_array = alignment_image(:,:,1); clear alignment_image;
        vasculature_confocal = rot90(imread('/Volumes/Data/2015-11-09-10/Imaging/confocal/C3-2015-11-09-10_tilestack_MIP.tif'));
        tubulin = rot90(imread('/Volumes/Data/2015-11-09-10/Imaging/confocal/C1-2015-11-09-10_tilestack_MIP.tif'));
        dapi = rot90(imread('/Volumes/Data/2015-11-09-10/Imaging/confocal/C2-2015-11-09-10_tilestack_MIP.tif'));
      
        pna = vasculature_confocal;
     
end


%% Load the image using Matlab's control point select tool. The confocal 
% image is the "moving" image, or the one that will rotate to conform to 
% the array image, the "fixed" image

if ~exist('movingPoints','var')
    cpselect(vasculature_confocal, vasculature_array);
    % After selecting the relevant pairs, File>Export Points to Workspace
end

%% Compute transform based on control points
% Specify the type of transformation and infer its parameters, using 
% fitgeotrans. fitgeotrans is a data-fitting function that determines the 
% transformation needed to bring the image into alignment, based on the 
% geometric relationship of the control points. fitgeotrans returns the 
% parameters in a geometric transformation object.
% Transformation type can be: 'nonreflectivesimilarity' | 'similarity' | 'affine' | 'projective'

% Piecewise linear transformation looked the best
tform_pwl = fitgeotrans(movingPoints,fixedPoints,'pwl');
registered_pwl = imwarp(vasculature_confocal, tform_pwl);
%% Display results
figure; imshow(registered_pwl); title('registered confocal vasculature image'); 
if strcmp(experiment_id, '2015-10-06-3')
    vasculature_array = padarray(vasculature_array,fliplr([3586 5923] - [3464 4418]),'pre');
elseif strcmp(experiment_id, '2015-11-09-3')
    registered_pwl = padarray(registered_pwl,[920 0],'pre');
    registered_pwl = registered_pwl(20:end,:); 
end
figure(100); imshow(vasculature_array); title('wide-field vasculature image over array');

%%
figure; imshow(registered_pwl);
hold on; h=imshow(vasculature_array); title('co-registered overlay using piecewise linear transformation');
hold off;

[M,N] = size(vasculature_array); 
block_size = 100; 
P = ceil(M / block_size); 
Q = ceil(N / block_size); 
alpha = checkerboard(block_size, ... 
    P, Q) > 0; 
alpha = alpha(1:M, 1:N); 
set(h, 'AlphaData', alpha);

%%
figure; 
r_im = uint8(zeros([size(registered_pwl) 3]));
r_im(:,:,1) = registered_pwl; 
g_im = uint8(zeros([size(vasculature_array) 3]));
row_shift = 0;
g_im(1:end-row_shift,:,2) = vasculature_array((row_shift+1):end,:); 
h1 = imshow(r_im);
hold on; h2=imshow(g_im); title('co-registered overlay using piecewise linear transformation');
hold off;

set(h1, 'AlphaData', 0.5);
set(h2, 'AlphaData', 0.33);
 %% % Pad the arrays so that the can be overlayed in RGB
% [morepadding] = size(vasculature_array) - size(registered_pwl); 
% array_padded = padarray(vasculature_array,[-1*morepadding(1) 0],'post'); 
% vasc_padded = padarray(registered_pwl,[0 morepadding(2)],'post'); 
% 
% rgb_merge = zeros([size(vasc_padded) 3]); 
% rgb_merge(:,:,1) = (array_padded - min(array_padded(:)))/(max(array_padded(:)) - min(array_padded(:)))*505; 
% rgb_merge(:,:,2) = (vasc_padded - min(vasc_padded(:)))/(max(vasc_padded(:))-min(vasc_padded(:)))*505; 
% figure; imshow(rgb_merge);

%% Plot the electrodes in proper image locations

if ~exist('newXYCoords','var')
    [xc,yc] = getElectrodeCoords512();
    yc = -yc;
    figure(100);
    [xx,yy] = ginput(4); % User clicks on the four corner electrodes.
    range_x = max(xc) - min(xc);
    range_xx = max(xx) - min(xx);
    yy_sort = sort(yy);
    yy = mean(reshape(yy_sort,2,[]));
    range_y = max(yc) - min(yc);
    range_yy = max(yy) - min(yy);
    
    scaleFactor = [range_yy/range_y range_xx/range_x];%mean(0.995*[range_yy/range_y range_xx/range_x]);
    offset = [repmat(xx(1),1,512);repmat(yy(1),1,512)];
    newXYCoords = [xc - min(xc); yc - min(yc)].*repmat(scaleFactor',1,512) + offset;
end
hold on;  scatter(newXYCoords(1,:),newXYCoords(2,:),30,[1 1 1], 'filled');
%% Warp the tubulin and DAPI images, Overlap with electrodes.

registered_tubulin = imwarp(tubulin,tform_pwl); 
registered_dapi = imwarp(dapi,tform_pwl); 
registered_vasc = imwarp(pna,tform_pwl); 
if strcmp(experiment_id, '2015-11-09-3')
    registered_tubulin = padarray(registered_tubulin,[920 0],'pre');
    registered_tubulin = registered_tubulin(20:end,:);
    registered_dapi = padarray(registered_dapi,[920 0],'pre');
    registered_dapi = registered_dapi(20:end,:);
    registered_vasc = padarray(registered_vasc,[920 0],'pre');
    registered_vasc = registered_vasc(20:end,:);
end

%% Plotting various ways
figure; imshow(registered_tubulin);
hold on;  scatter(newXYCoords(1,:),newXYCoords(2,:),30,[1 1 0], 'filled');
figure; imshow(registered_dapi);
hold on;  scatter(newXYCoords(1,:),newXYCoords(2,:),30,[1 1 0], 'filled');

figure; imshow(registered_pwl);
hold on;  scatter(newXYCoords(1,:),newXYCoords(2,:),30,[1 1 0], 'filled');
for e = 1:512; 
    text(newXYCoords(1,e),newXYCoords(2,e)-50,num2str(e),'HorizontalAlignment','center','Color',[0 1 1]); 
end

rgb_merge = cat(3,registered_pwl,registered_tubulin,registered_dapi); 
figure; imshow(rgb_merge);
hold on;  scatter(newXYCoords(1,:),newXYCoords(2,:),30,[1 1 1], 'filled');
for e = 1:512; 
    text(newXYCoords(1,e),newXYCoords(2,e)-50,num2str(e),'HorizontalAlignment','center','Color',[1 1 1]); 
end
bg_only = rgb_merge; 
bg_only(:,:,1) = 0; 
figure; imshow(bg_only); 
rg_only = rgb_merge; 
rg_only(:,:,3) = 0; 
figure; imshow(rg_only); 
g_only = rgb_merge; 
g_only(:,:,[1 3]) = 0; 
figure; imshow(g_only)



%% Other tries
% mytform = fitgeotrans(movingPoints1, fixedPoints1, 'projective');
% tform_pwl = fitgeotrans(movingPoints1, fixedPoints1,'pwl');
% 
% mytform_a = fitgeotrans(movingPoints1, fixedPoints1, 'affine');
% mytform_s = fitgeotrans(movingPoints1, fixedPoints1, 'similarity');
% 
% mytform_lwm = fitgeotrans(movingPoints, fixedPoints,'lwm',6);
% registered = imwarp(vasculature_array, mytform);
% registered_pwl = imwarp(vasculature_array, tform_pwl); % This one is better
% registered_a = imwarp(vasculature_array, mytform_a); %worse than pwl
% registered_s = imwarp(vasculature_array, mytform_s);
% registered_lwm = imwarp(vasculature_array, mytform_lwm);
% figure; imshow(registered); 
% figure; imshow(registered_pwl); 
% mask = zeros(size(registered_pwl)); 