%% Script to analyze axon bundle images. 
% Need to have loaded the "registered_tubulin" image from
% registerArrayToConfocal.m, and newXYCoords

% OR %
%Try loading the images edited in photoshop.

im1 = imread('/Volumes/Analysis/2015-10-06-6/image analysis/photoshop_warp/array.jpg'); 
im2 = imread('/Volumes/Analysis/2015-10-06-6/image analysis/photoshop_warp/tubulin_warped.jpg'); 
im3 = imread('/Volumes/Analysis/2015-10-06-6/image analysis/photoshop_warp/vasculature_warped.jpg'); 
array = im1(817:5288,393:8968,1);
for ii = 1:size(array,1)
    if all(array(ii,:) ==255)
        disp(['discard row ' num2str(ii)]); 
    end
end
vasc = im3(:,:,1);
B = imresize(array,size(vasc)); 
rg = uint8(zeros([size(B) 3])); 
rg(:,:,1)= B; 
rg(:,:,2) = vasc*0.7; 
figure; imshow(rg); 
figure; imshow(im3(:,:,1)); 
maskOfElectrodes =uint8(zeros(size(B))); 
maskOfElectrodes(B < 5) = 255; 
figure; imshow(maskOfElectrodes)
rg(:,:,3) = maskOfElectrodes; 
figure; imshow(rg); 

figure; imshow(im2); 
tubulin = uint8(zeros(size(im2))); 
tubulin(:,:,2) = im2(:,:,1) + maskOfElectrodes; 
tubulin(:,:,1 ) = maskOfElectrodes; 
tubulin(:,:,3) = maskOfElectrodes; 
figure; imshow(tubulin); 
hsize = 60% 120; % 241=60 um separation of electrodes.
sigma = 0.46;
h = fspecial('log',hsize,sigma);
filtered = imfilter(im2(:,:,1),h);
figure; imshow(filtered); title(sprintf('laplacian of a gaussian hsize %0.0f sigma %0.2f',hsize,sigma)) ;


[xc,yc] = getElectrodeCoords512();
if exist('newXYCoords','var') == 0
    yc = -yc;
    [xx,yy] = ginput(4); % User clicks on the four corner electrodes.
    range_x = max(xc) - min(xc);
    range_xx = max(xx) - min(xx);
    yy_sort = sort(yy);
    yy = mean(reshape(yy_sort,2,[]));
    range_y = max(yc) - min(yc);
    range_yy = max(yy) - min(yy);
    
    scaleFactor = 0.9995*[range_yy/range_y range_xx/range_x];%mean(0.995*[range_yy/range_y range_xx/range_x]);
    offset = [repmat(xx(1),1,512);repmat(yy(1),1,512)];
    newXYCoords = [xc - min(xc); yc - min(yc)].*repmat(scaleFactor',1,512) + offset;
end
hold on;  scatter(newXYCoords(1,:),newXYCoords(2,:),30,[1 1 0], 'filled');

load('/Users/grosberg/Lab/matlab/code/projects/electrical_stim/resources/array_matrix_id510.mat')
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/axonBundleThresholds_2015_10_06_6.mat'); 
idx = find(bundleThresholds_2015_10_06_6);
idx0 = find(bundleThresholds_2015_10_06_6 == 0);

% overlay with tubulin stain 
hold on; scatter(newXYCoords(1,idx),newXYCoords(2,idx),100*bundleThresholds_2015_10_06_6(idx),'r','filled'); axis image; axis off; colormap(jet);
figure; imshow(filtered); 
hold on; scatter(newXYCoords(1,idx),newXYCoords(2,idx),100*bundleThresholds_2015_10_06_6(idx),'r','filled'); axis image; axis off; 
title('Bundle thresholds by electrodes - larger dot = higher threshold')

% Decide which electrodes to use (discord those near tears in the retina, for example). 
yheight = 64;

nxc = sort(unique(newXYCoords(1,:))); 
nyc = sort(unique(newXYCoords(2,:))); 
[xc,yc] = getElectrodeCoords512();

figure; scatter(xc(idx),yc(idx),100,bundleThresholds_2015_10_06_6(idx),'filled'); axis image; axis off; colormap(jet); colorbar; 
hold on; scatter(xc(idx0),yc(idx0),100,'k'); 
title('Bundle threshold values'); 
figure(10); 
for row = 4:15
    cla; 
    elecs = array_matrix_id510(row,:); 
    elecs_subset = elecs(5:22);  
   figure(5);
   hold on; 
   for e = 1:length(elecs)
       text(newXYCoords(1,elecs(e)),newXYCoords(2,elecs(e)),num2str(elecs(e)),...
           'HorizontalAlignment','center','FontSize',18,'Color','g');
   end
    figure(10); 
    [AX,H1,H2] =plotyy(nxc(1):nxc(end),mean(filtered((nyc(row)-yheight/2):(nyc(row)+yheight/2),nxc(1):nxc(end)),1),newXYCoords(1,elecs),bundleThresholds_2015_10_06_6(elecs));
    set(H2(1),'Marker','d','MarkerSize',14,'LineStyle','none')
    title(['row ' num2str(row)]); 
    pause; 
end
%%
fluo_vals = []; 
bundle_ts = [];
figure; 
for row = 4:12
    elecs = array_matrix_id510(row,:);  
    elecs = elecs(5:22);  
    image_pixels = nxc(1):nxc(end);
    [Lia,Locb] = ismember(round(newXYCoords(1,elecs)),round(image_pixels));
    image_selection = mean(filtered((nyc(row)-yheight/2):(nyc(row)+yheight/2),nxc(1):nxc(end)),1);
    if all(Lia)
        vals = image_selection(Locb);
    else
        keyboard; 
    end
    hold on; 
    scatter(vals,bundleThresholds_2015_10_06_6(elecs),20,'r','filled');
    %     title(['row ' num2str(row)]);
    fluo_vals = [fluo_vals vals];
    bundle_ts= [bundle_ts bundleThresholds_2015_10_06_6(elecs)];
end

[r,p] = corrcoef(fluo_vals,bundle_ts);
title(sprintf('Correlation is %0.4f with p value of %0.3f',r(2),p(2))); 
xlabel('fluoresence intensity value'); 
ylabel('bundle activation threshold'); 

%% Create an axon bundle mask. 
maskOfBundles =uint8(zeros(size(filtered))); 
maskOfBundles(filtered >120) = 255; 
figure; imshow(maskOfBundles);

%%
radius = 120; %pixels, 120 is equal to 30 microns. 
[row_idx,col_idx]= meshgrid(1:size(maskOfBundles,2),1:size(maskOfBundles,1)); 

figure;
fluo_vals = zeros(1,162);
bundle_ts = zeros(1,162); 
patternNos = zeros(1,162); 
ii = 1; 
for row = 4:12
    elecs = array_matrix_id510(row,:);
    elecs = elecs(5:22);
    hold on;
    for e = 1:length(elecs)
        %         text(newXYCoords(1,elecs(e)),newXYCoords(2,elecs(e)),num2str(elecs(e)),...
        %             'HorizontalAlignment','center','FontSize',18,'Color','g');
        % Find pixels within a given radius.
        C = sqrt((row_idx-newXYCoords(1,elecs(e))).^2 + ...
            (col_idx-newXYCoords(2,elecs(e))).^2)<=radius;
%         hold on; scatter(mean(maskOfBundles(C)),bundleThresholds_2015_10_06_6(elecs(e)),30,'r','filled');
        hold on; scatter(mean(filtered(C)),bundleThresholds_2015_10_06_6(elecs(e)),30,'r','filled');
        patternNos(ii) = elecs(e); 
        fluo_vals(ii) = mean(filtered(C)); 
        bundle_ts(ii) = bundleThresholds_2015_10_06_6(elecs(e)); 
        ii = ii + 1; 
        disp(['Finished with electrode ' num2str(elecs(e))]);
    end
end
xlabel(sprintf('Mean intensity in a %0.0f um radius',radius/120*30)); 
ylabel('axon bundle activation threshold (\muV)'); 

[r,p] = corrcoef(fluo_vals,bundle_ts);
title(sprintf('Correlation is %0.4f with p value of %0.3f',r(2),p(2))); 
figure; scatter(fliplr(fluo_vals),bundle_ts,30,'b','filled');
[r,p] = corrcoef(fliplr(fluo_vals),bundle_ts);
title(sprintf('Fluo vals flipped \ncorrelation is %0.4f with p value of %0.3f',r(2),p(2)));
xlabel(sprintf('Mean intensity in a %0.0f um radius',radius/120*30)); 
ylabel('axon bundle activation threshold (\muA)'); 

%%
radius = [5 10 15 20 25 30]; 
corr_val = [-0.2156 -0.2575 -0.3272 -0.3083 -0.2150 -0.1463]; 
figure; scatter(radius, abs(corr_val),100,'k','filled'); 
xlabel('Radius (um) over which fluorescence averaged'); 
ylabel('Absolute value of correlation coeffienct'); 
%% Find maximum recorded signal for the bundles. 
pathName = '/Volumes/Analysis/2015-10-06-6/data001/';
display = 0; 
exclude = 0; 
bundleMeans = getBundleVoltagesAStar(pathName, patternNos, display, exclude); 

figure; 
for n = 1:size(bundleMeans,3)
    idx = find(bundleMeans(:,1,n));
    hold on; plot(abs(bundleMeans(idx,2,n)),abs(bundleMeans(idx,1,n))); 
end
avg_rec_sig = max(abs(squeeze(bundleMeans(:,1,:))),[],1);
%% make some plots. 
figure; scatter(avg_rec_sig,bundle_ts,30,'b','filled');
xlabel('mean recorded voltage over bundle at max stim amp (uV)'); 
ylabel('bundle activation threshold (uA)'); 
[r,p] = corrcoef(avg_rec_sig,bundle_ts);
title(sprintf('correlation is %0.4f with p value of %0.3f',r(2),p(2)));

figure; scatter(fliplr(fluo_vals),avg_rec_sig,30,'b','filled');
[r,p] = corrcoef(fliplr(fluo_vals),avg_rec_sig);
title(sprintf('Fluo vals flipped \ncorrelation is %0.4f with p value of %0.3f',r(2),p(2)));
xlabel(sprintf('Mean intensity in a %0.0f um radius',radius/120*30)); 
ylabel('mean recorded voltage over bundle at max stim amp (uV)'); 

figure; scatter((fluo_vals),avg_rec_sig,30,'r','filled');
[r,p] = corrcoef((fluo_vals),avg_rec_sig);
title(sprintf('correlation is %0.4f with p value of %0.3f',r(2),p(2)));
xlabel(sprintf('Mean intensity in a %0.0f um radius',radius/120*30)); 
ylabel('mean recorded voltage over bundle at max stim amp (uV)'); 

figure; imshow(filtered); 
hold on; scatter(newXYCoords(1,patternNos),newXYCoords(2,patternNos),...
    2*avg_rec_sig,'r','filled'); axis image; axis off; 
title('Average recorded signal from bundles activated by electrode'); 

figure; imshow(filtered); 
hold on; scatter(newXYCoords(1,patternNos),newXYCoords(2,patternNos),...
    2*fluo_vals+0.1,'g','filled'); axis image; axis off; 
title(sprintf('Mean intensity in a %0.0f um radius',radius/120*30));
%% Plotting intensity profiles across electrodes. 
if 0
%     smaller = registered_tubulin(750:3200,1:3200);
    hsize = 120; % 241=60 um separation of electrodes. 
    sigma = 0.46; 
    h = fspecial('log',hsize,sigma);
    filtered = imfilter(registered_tubulin,h); 
    figure; imshow(filtered); title(sprintf('laplacian of a gaussian hsize %0.0f sigma %0.2f',hsize,sigma)) ;
    figure; imshow(registered_tubulin);
    hold on;  scatter(newXYCoords(1,:),newXYCoords(2,:),30,[1 1 0], 'filled');
end
[x, y] = ginput(4);
yheight = 64;

% figure; imshow(registered_tubulin(1167:1167+yheight,1400:8978)); 
% figure; plot(mean(registered_tubulin(1167:1167+yheight,1400:8978),1));
nxc = sort(unique(newXYCoords(1,:))); 
nyc = sort(unique(newXYCoords(2,:))); 


[xc,yc] = getElectrodeCoords512();
load('/Users/grosberg/Lab/matlab/code/projects/electrical_stim/resources/array_matrix_id510.mat')
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/axonBundleThresholds_2015_10_06_6.mat'); 
idx = find(bundleThresholds_2015_10_06_6);
idx0 = find(bundleThresholds_2015_10_06_6 == 0);
figure; scatter(xc(idx),yc(idx),100,bundleThresholds_2015_10_06_6(idx),'filled'); axis image; axis off; colormap(jet); colorbar; 
hold on; scatter(xc(idx0),yc(idx0),100,'k'); 
% overlay with tubulin stain 
hold on; scatter(newXYCoords(1,idx),newXYCoords(2,idx),100*bundleThresholds_2015_10_06_6(idx),'g','filled'); axis image; axis off; colormap(jet);
figure; 
for row = 4:16
    cla; 
    elecs = array_matrix_id510(row,:); 
    [AX,H1,H2] =plotyy(nxc(1):nxc(end),mean(filtered((nyc(row)-yheight/2):(nyc(row)+yheight/2),nxc(1):nxc(end)),1),newXYCoords(1,elecs),bundleThresholds_2015_10_06_6(elecs));
    set(H2(1),'Marker','d','MarkerSize',14,'LineStyle','none')
    title(['row ' num2str(row)]); 
    pause; 
end

figure; 
colors = colorcube(16); 
for row = 4:16
    elecs = array_matrix_id510(row,:);   
    image_pixels = nxc(1):nxc(end);
    [Lia,Locb] = ismember(round(newXYCoords(1,elecs)),round(image_pixels));
    image_selection = mean(filtered((nyc(row)-yheight/2):(nyc(row)+yheight/2),nxc(1):nxc(end)),1);
    if all(Lia)
        vals = image_selection(Locb);
    else
        keyboard; 
    end
    hold on; 
    scatter(vals,bundleThresholds_2015_10_06_6(elecs),20,'r','filled');
%     title(['row ' num2str(row)]); 
end
xlabel('fluorescence level'); 
ylabel('axon bundle activation threshold'); 
fimage = imread('/Users/grosberg/Desktop/filteredtubulin2.jpg'); 

%% Analyze 2015-11-09-3
im1 = imread('/Volumes/Analysis/2015-11-09-3/image analysis/registered_tubulin_crop.tiff');
load('/Volumes/Analysis/2015-11-09-3/image analysis/registered_electrode_coordinates_crop.mat'); 
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/AlgorithmOutputs/axonBundleThresholdsAlgorithm_2015_11_09_3.mat')
Safe1193 = axonBundleThresholdsAlgorithm_2015_11_09_3;
    
hsize = 120; % 241=60 um separation of electrodes.
sigma = 0.42;
h = fspecial('log',hsize,sigma);
filtered = imfilter(im1,h);
figure; imshow(filtered); title(sprintf('laplacian of a gaussian hsize %0.0f sigma %0.2f',hsize,sigma)) ;


idx = find(Safe1193~=4); 
elecs = 1:512; 
%%
radii = [5 10 15 20 30 40 60 120 240];
cor_coeffs = zeros(size(radii));
pvals = zeros(size(radii));
cor_coeffs2 = zeros(size(radii));
pvals2 = zeros(size(radii));
for r = 1:length(radii)
    radius = radii(r); %pixels, 120 is equal to 30 microns.
    [row_idx,col_idx]= meshgrid(1:size(filtered,2),1:size(filtered,1));
    fluo_vals = zeros(512,1);
    for e = 1:length(elecs)
        % Find pixels within a given radius.
        C = sqrt((row_idx-newXYCoords_crop(1,elecs(e))).^2 + ...
            (col_idx-newXYCoords_crop(2,elecs(e))).^2)<=radius;
        fluo_vals(e) = mean(filtered(C));
        disp(['Finished with electrode ' num2str(elecs(e))]);
    end
    [cc,p]=corrcoef(fluo_vals(idx),Safe1193(idx));
    figure; scatter(fluo_vals(idx),Safe1193(idx),10,'k','filled');
    xlabel(sprintf('Mean intensity in a %0.0f um radius',radius/120*30));
    ylabel('axon bundle activation threshold (\muV)');
    title(sprintf('Correlation coefficient is %0.4f, p=%0.2f',cc(2),p(2)));
    cor_coeffs(r) = cc(2);
    pvals(r) = p(2);
    [cc,p]=corrcoef(fluo_vals(analyzedElecs),safeDiff(analyzedElecs));
    figure; scatter(fluo_vals(analyzedElecs),safeDiff(analyzedElecs),10,'k','filled');
    xlabel(sprintf('Mean intensity in a %0.0f um radius',radius/120*30));
    ylabel('somatic - axon bundle activation threshold (\muV)');
    title(sprintf('Correlation coefficient is %0.4f, p=%0.2f',cc(2),p(2)));
    cor_coeffs2(r) = cc(2);
    pvals2(r) = p(2);
end

figure; plot(radii/120*30,cor_coeffs,'-x');
ylabel('correlation coefficient');
title('correlation between bundle threshold & fluorescence intensity'); 
xlabel('radius around electrode (um) for fluorescence averaging');