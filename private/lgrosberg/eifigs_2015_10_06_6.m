
playMovie512arrayAfterStimPattern_dots('/Volumes/Analysis/2015-10-06-6/data001/', 308, 'movieIndex', 25);
playMovie512arrayAfterStimPattern_dots('/Volumes/Analysis/2015-10-06-6/data001/', 308, 'movieIndex', 27);

[~, amps0] = generateEiFromStimPattern('/Volumes/Analysis/2015-10-06-6/data001/',...
    308,'movieNo',194,'plotElecWaveforms',[73 69 59 55 52 282 278 269],'showElecNums',false);
[~, amps1] = generateEiFromStimPattern('/Volumes/Analysis/2015-10-06-6/data001/', 308,'movieNo',202);

[~, amps0] = generateEiFromStimPattern('/Volumes/Analysis/2015-10-06-6/data001/',...
    308,'movieNo',202,'plotElecWaveforms',[73 69 59 55 52 282 278 269],'showElecNums',false);
% There is a slight axonal activation at the lower current amplitude. 

%%
playMovie512arrayAfterStimPattern_dots('/Volumes/Analysis/2015-10-06-6/data001/', 316, 'movieIndex', 25);
playMovie512arrayAfterStimPattern_dots('/Volumes/Analysis/2015-10-06-6/data001/', 308, 'movieIndex', 27);
playMovie512arrayAfterStimPattern_dots('/Volumes/Analysis/2015-10-06-6/data001/', 86, 'movieNo', 299);


playMovie512arrayAfterStimPattern_dots('/Volumes/Analysis/2015-10-06-6/data001/', 90, 'movieNo', 197);
[rawdata, amps0] = generateEiFromStimPattern('/Volumes/Analysis/2015-10-06-6/data001/',...
    90,'movieNo',197);
eiContour_wLinFit(amps0);
alltrials = 1:25; 
spiketrials = [2 3 4 6 10 11 14 15 16 17 18 19 21 23];
nospiketrials = setxor(alltrials,spiketrials);
figure; 
toplot = squeeze(rawdata(:,93,:))';
subM = repmat(mean(toplot,2),1,size(toplot,2)); 
data = toplot-subM;
figure; 
plot(data(:,spiketrials),'r'); 
hold on; plot(data(:,nospiketrials),'k');

spikedata = squeeze(mean(rawdata(spiketrials,:,:),1)); 
subdata = squeeze(mean(rawdata(nospiketrials,:,:),1));
data = spikedata - subdata; 
amps = max(data,[],2) - min(data,[],2);
eiContour_wLinFit(amps)

playMovie512arrayAfterStimPattern_dots('/Volumes/Analysis/2015-10-06-6/data001/', 90, 'movieNo', 221);
[rawdata, amps0] = generateEiFromStimPattern('/Volumes/Analysis/2015-10-06-6/data001/',...
    90,'movieNo',221);


[rawdata, amps0] = generateEiFromStimPattern('/Volumes/Analysis/2015-10-06-6/data001/',...
    90,'movieNo',197,'plotElecWaveforms',[238 230 231 232 259 264 276 105 93 86 83 79 76 72 64]);
[rawdata, amps0] = generateEiFromStimPattern('/Volumes/Analysis/2015-10-06-6/data001/',...
    90,'movieNo',221,'plotElecWaveforms',[238 230 231 232 259 264 276 105 93 86 83 79 76 72 64]);

% Contour plot of n1338
[eiM,neuronIdList] = convertEiFileToMatrix('/Volumes/Analysis/2015-10-06-6/data002/data002.ei');
ei = eiM(:,:,find(neuronIdList==1338)); 
eiamps = max(ei,[],2) - min(ei,[],2); 

eiContour_wLinFit(eiamps)
patternNos = find(positions(:,1)< max(positions(:,1))/2 & ...
    positions(:,1)> -max(positions(:,1))/2 &...
    positions(:,2)< max(positions(:,2))/2 & ...
    positions(:,2)> -max(positions(:,2))/2);
figure; 
scatter(positions(:,1),positions(:,2)); 
hold on; scatter(positions(patternNos,1),positions(patternNos,2),'k');
figure;
for p = 1:length(patternNos)
    bundleMeans = getBundleVoltagesAStar('/Volumes/Analysis/2015-10-06-6/data001/', ...
        patternNos(p), 0, 0);
    hold on;
    plot(abs(bundleMeans(1:38,2)),abs(bundleMeans(1:38,1)));
    disp(['finished with pattern ' num2str(p) ' of ' num2str(length(patternNos))]); 
end
xlabel('stimulation amplitude (uA)'); ylabel('mean voltage recorded along bundle'); 

figure; 
scatter(abs(bundleMeans(1:38,2)),abs(bundleMeans(1:38,1)),100,'k','filled');
xlabel('stimulation amplitude (uA)'); ylabel('mean voltage recorded along bundle'); 
stimamps = abs(bundleMeans(1:38,2));
data = abs(bundleMeans(1:38,1)); 
cummean = cumsum(data)./(1:length(data))'; 
figure; plot(stimamps,cummean,'-d'); title('cumulative mean'); 
figure; plot(diff(cummean),'-d'); title('change in cumulative mean'); 

% deviation of the next point from the previous
sq_err = zeros(length(data),1); 
for j = 2:length(data); 
    sq_err(j) = (data(j) - cummean(j-1))^2; 
end
idx = find(sq_err>50,1);
figure; plot(stimamps,sq_err,'-d'); hold on; 
plot(stimamps(idx),sq_err(idx),'ro'); 

% Use to calculate the bundle activation threshold. 
patternNos = 1:512; 
err_thresh = 50; % may want to optimize this value. 
bundleThresh = zeros(length(patternNos),1); 
for p = 56:length(patternNos)
    bundleMeans = getBundleVoltagesAStar('/Volumes/Analysis/2015-10-06-6/data001/', ...
        patternNos(p), 0, 0);
    data = abs(bundleMeans(1:38,1)); 
    cummean = cumsum(data)./(1:length(data))';
    % deviation of the next point from the previous
    sq_err = zeros(length(data),1);
    for j = 2:length(data);
        sq_err(j) = (data(j) - cummean(j-1))^2;
    end
    idx = find(sq_err > err_thresh,1);
    if isempty(idx) % no bundle activation found
        bundleThresh(p) = 4;
    else
        bundleThresh(p) = abs(bundleMeans(idx,2));
    end
    disp(['finished with pattern ' num2str(p) ' of ' num2str(length(patternNos))]); 
end
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/hand_sorted_bundle_thresholds/axonBundleThresholds_2015_10_06_6.mat')
figure; scatter(bundleThresh,bundleThresholds_2015_10_06_6,10,'k'); 
xlabel('A* algorithm bundle threshold'); 
ylabel('human bundle threshold'); title('2016-10-06-6');
hold on; line(0:4,0:4); 
% translated histogram
points = [bundleThresholds_2015_10_06_6' ; bundleThresh'];
dot(bundleThresholds_2015_10_06_6,10,bundleThresh)
linecoords = ones(size(points)); 
linecoords(2,:) = -1*linecoords(2,:);
proj = dot(points, linecoords)./dot(linecoords,linecoords)*sqrt(2); 
figure; plot(proj);
x =[-4:0.25:4]; 
figure; subplot(1,2,1); scatter(bundleThresh,bundleThresholds_2015_10_06_6,10,'k'); 
xlabel('A* algorithm bundle threshold'); ylabel('human bundle threshold'); title('2016-10-06-6');
hold on; line(0:4,0:4); 
subplot(1,2,2); hist(proj,x);xlim([-4.5 4.5]); xlabel('human - algorithm'); ylabel('electrodes'); title('2016-10-06-6');
% figure; hist(proj,x); ylim([0 28]); xlim([-3 3]);

% Plot the current spread of the signal over the array.
[xc,yc] = getElectrodeCoords512(); 
[rawData, amplitudes] = generateEiFromStimPattern('/Volumes/Analysis/2015-10-06-6/data001/',...
    1,'suppressPlots',true); 

firstArtifact = mean(rawData(:,:,:,1),1);
subtractionMatrix = repmat(firstArtifact,[size(rawData,1) 1]);
subtractionMatrix =repmat(subtractionMatrix,[1 1 1 size(rawData,4)]); 
modData = rawData - subtractionMatrix;
meanData = squeeze(mean(modData,1)); 
allAmps = squeeze(max(meanData,[],2) - min(meanData,[],2)); 
idx = find(yc==210); 
[~,ii]=sort(xc(idx));
figure; plot(allAmps(idx,38)); hold on; plot(amplitudes(idx),'-o'); 
figure; plot(allAmps(idx,:)); 

%
xq = min(xc):30:max(xc); 
yq = min(yc):30:max(yc); 
vq = zeros(size(yq,2),size(xq,2),size(allAmps,2)); 
for a = 1:size(allAmps,2)
    vq(:,:,a) = griddata(xc,yc,allAmps(:,a),xq,yq');
end
figure; imagesc(xq,yq,vq); axis xy; axis image; 
F = scatteredInterpolant(xc',yc',allAmps(:,38)); 
%% Show movies of stimulation overlayed with the tubulin image

% Lower amplitude activation, shows the activation of one cell. 
playMovie512arrayAfterStimPattern_dots('/Volumes/Analysis/2015-10-06-6/data001/', 90, 'movieNo', 197);

% Axon bundle activation at a higher stimulus amplitude. 
playMovie512arrayAfterStimPattern_dots('/Volumes/Analysis/2015-10-06-6/data001/', 90, 'movieNo', 221);

% Load tubulin image. 
axon_img = imread(['/Volumes/Analysis/2015-10-06-6/image analysis/'...
    'photoshop_warp/tubulin_warped.jpg']); 
hsize = 120; % 241=60 um separation of electrodes.
sigma = 0.46;
h = fspecial('log',hsize,sigma);
filtered = imfilter(im2(:,:,1),h);

% Load electrode coordinates. 
load(['/Volumes/Analysis/2015-10-06-6/image analysis/photoshop_warp/'...
    'electrode_coords.mat']);
pathToAnalysisData = '/Volumes/Analysis/2015-10-06-6/data001/'; 
patternNo = 90; 
movieNos = findMovieNos(pathToAnalysisData, patternNo); 
playMovie512arrayAfterStimPattern_imageOverlay(pathToAnalysisData, patternNo,...
    filtered, newXYCoords, 'movieNo', 221,'saveMovie',true); 

playMovie512arrayAfterStimPattern_imageOverlay(pathToAnalysisData, patternNo,...
    filtered, newXYCoords, 'movieNo', 197,'saveMovie',true); 

playMovie512arrayAfterStimPattern_imageOverlay(pathToAnalysisData, patternNo,...
    filtered, newXYCoords, 'movieNo', 309,'saveMovie',true); 

%% Validate the polynomial fits
[eiM,neuronIdList] = convertEiFileToMatrix('/Volumes/Analysis/2015-10-06-6/data002/data002.ei');
[xc,yc] = getElectrodeCoords512();
figure; 
for n = 1: length(neuronIdList)
    cla;
    for i = 1:length(xc)
        plot(xc(i), yc(i), 'ok','MarkerSize', 2, 'MarkerFaceColor', 'k')
    end
    eiAmps = max(eiM(:,:,n),[],2) - min(eiM(:,:,n),[],2);
    scatter(xc,yc,10*eiAmps+0.1,'b','filled'); 
    [xf, yf] = weighted_axon_poly_reg(eiAmps, 'ei_thresh', 4);
    hold on; plot(xf,yf,'-','Color','green','LineWidth',2);
    axis image; axis off;   
    title(sprintf('neuron id %0.0f',neuronIdList(n))); 
    pause; 
end

% Load tubulin image. 
axon_img = imread(['/Volumes/Analysis/2015-10-06-6/image analysis/'...
    'photoshop_warp/tubulin_warped.jpg']); 
hsize = 120; % 241=60 um separation of electrodes.
sigma = 0.46;
h = fspecial('log',hsize,sigma);
filtered = imfilter(axon_img(:,:,1),h);

% Load electrode coordinates. 
load(['/Volumes/Analysis/2015-10-06-6/image analysis/photoshop_warp/'...
    'electrode_coords.mat']);
scalef = mean((max(newXYCoords')-min(newXYCoords'))./[max(xc)-min(xc),max(yc)-min(yc)]);  
offsetx = mean(newXYCoords(1,:) - scalef*xc); 
offsety = mean(newXYCoords(2,:) + scalef*yc);

%% And now with the image overlay
    [movfile,movpath] = uiputfile('*.avi','Save Movie As');
    movieFileName = [movpath movfile];
    writerObj = VideoWriter(movieFileName);
    writerObj.FrameRate = 3;
    open(writerObj);
f = figure; set(f,'Position',[100 465 845 445]);
set(f,'Color','white');
imshow(filtered); axis image; 
% % Show electrode positions. 
% hold on; scatter(newXYCoords(1,:),newXYCoords(2,:),30,[1 1 0],'filled'); 

for n = 1: length(neuronIdList)
    % Get EI of current neuron 
    eiAmps = max(eiM(:,:,n),[],2) - min(eiM(:,:,n),[],2);
    % Fit a polynomial regression
    [xf, yf] = weighted_axon_poly_reg(eiAmps, 'ei_thresh', 4);
    % Scatter plot of the EI
    hold on;  hei = scatter(newXYCoords(1,:),newXYCoords(2,:),20*eiAmps+0.1,'y','filled'); 
    % Overlay the polynomial fits
    hold on; h = plot(xf*scalef+offsetx,-yf*scalef+offsety,'-','Color','green','LineWidth',3);
    title(sprintf('neuron id %0.0f',neuronIdList(n))); 
    pause(0.01);  
    M = getframe(f);writeVideo(writerObj,M);
    delete(hei);          
    pause(0.01);  
    M = getframe(f);writeVideo(writerObj,M);
    delete(h);
    M = getframe(f);writeVideo(writerObj,M);
end

close(writerObj);