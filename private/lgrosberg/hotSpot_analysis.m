%% Find the percentage signal per electrode. 
eipath = '/Volumes/Analysis/2015-11-09-3/data000/data000.ei';
[eiM,nIdList] = convertEiFileToMatrix(eipath);
valsAllCells = zeros(size(eiM,1),size(eiM,3)); 

[xc,yc] = getElectrodeCoords512(); 
for c = 1:size(eiM,3)
    [mins,idx] = sort(min(eiM(:,:,c),[],2)); 
    scaledVals = mins/mins(1); 
    % Penalize edge cases.
    if abs(xc(idx(1)))>=915 || abs(yc(idx(1))) == max(yc)
        scaledVals(2) = 1; 
    end
    valsAllCells(:,c) = scaledVals; 
end
%%
secondaryElecs = 2:8; 
secondaryVals = sum(valsAllCells(secondaryElecs,:),1); 
[V,I] = sort(secondaryVals); 
sortedHotSpotIdx = 1./V; 
% Look at the first 10 neurons to see their 'hot spots'
nIdList(I(1:10));
eiAmps = squeeze(max(eiM,[],2) - min(eiM,[],2)); 

% figure; scatter(xc,yc,10,'k','filled'); axis image; axis off; 
fPath = '/Volumes/Analysis/2015-11-09-3/data001-data002-autosort/';

% Calculate the thresholds of all the cells ranked by hot-spottiness
allthresholds= zeros(1,size(eiAmps,2)); 
allmaxelecs = zeros(1,size(eiAmps,2)); 
% For display only
for ii = 1:4
    idxx = find(eiAmps(:,I(ii))>1);
    figure; subplot(1,2,1); 
    hold on; scatter(xc(idxx),yc(idxx),eiAmps(idxx,I(ii)),'filled'); axis image; axis off; 
    [axon_x, axon_y, ~, soma_x, soma_y] = weighted_axon_poly_reg(eiAmps(:,I(ii)));
    plot(axon_x, axon_y, 'Color', 'k','LineWidth',2);
    scatter(soma_x,soma_y,max(eiAmps(:,I(ii)))+0.1,'k','filled');  
    title(['n' num2str( nIdList(I(ii))) ]); 
    subplot(1,2,2); 
     [resort,~]  = sort(min(eiM(:,:,I(ii)),[],2)); 
    bar(abs(resort(1:10))); xlabel('electrode order of signal strength'); 
    ylabel('measured voltage'); title(['hot-spot value is ' num2str(sortedHotSpotIdx(ii))]);
end  
%%
for ii = 1:size(eiAmps,2)
    idxx = find(eiAmps(:,I(ii))>1);
    [resort,idx]  = sort(min(eiM(:,:,I(ii)),[],2)); 
    fname = ['elecRespAuto_n' num2str( nIdList(I(ii))) '_p' num2str(idx(1)) '.mat']; 
    try
        temp = load(fullfile(fPath,fname));
        [threshold,~, ~]  = getThresholdFromAuto(temp.elecRespAuto);
        if threshold>5; threshold = 5; end
        if threshold<0; threshold = 5; end
        allthresholds(ii) = threshold;
        allmaxelecs(ii) = idx(1);
    catch
        disp([fname ' does not exist; signal of this neuron is ' num2str(resort(1))]);
        allthresholds(ii) = nan;
        allmaxelecs(ii) = nan;
    end
end

load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/data_files/Safe1193.mat');
Safe1193(find(Safe1193==0)) = nan;
figure; plot(sortedHotSpotIdx,allthresholds,'kd'); title('Somatic activation thresholds'); 
ylabel('stim amp (uA)'); xlabel('Hot value 1/sum(scaled electrodes 2:5 (higher -> more hot-spotty)');
figure; hold on; 
for ii = size(allthresholds,2):-1:1
    if isfinite(allthresholds(ii))
        ratio(ii) = allthresholds(ii)/(Safe1193(allmaxelecs(ii))/1.1); 
        scatter(allthresholds(ii),Safe1193(allmaxelecs(ii))/1.1,50,sortedHotSpotIdx(ii),'filled')
    end
end
line(0:5,0:5,'Color','k');
xlabel('Somatic 50% activation threshold at electrode with max signal'); 
ylabel('Bundle safe zone,uA  (largest current with no bundle activation)'); 
c = colorbar; caxis([0 2]); colormap jet; 
c.Label.String = 'hot-spotty value 1/sum(electrodes 2:5'; 
figure; hold on;
for ii = 1:size(allthresholds,2)
    if isfinite(allthresholds(ii))
        ratio(ii) = allthresholds(ii)/(Safe1193(allmaxelecs(ii))/1.1); 
        if isfinite(ratio(ii))
            if ratio(ii) > 1
                a=area([ii-1 ii],ratio(ii)*[1 1]);
                a.FaceColor = [1 0 0];
                a.EdgeColor = [1 0 0];
            elseif ratio(ii)<1
                a=area([ii-1 ii],ratio(ii)*[1 1]);
                a.FaceColor = [0.5 0.9 0.5];
                a.EdgeColor = [0.5 0.9 0.5];
            end
        end
    end
end
xlim([1 size(allthresholds,2)]); 
xlabel('cell hot-spottiness ranking, 1 is highest')
ylabel('ratio: cell threshold : axon bundle threshold')
 
line(get(gca,'XLim'),[1 1],'Color','k','LineStyle','--')
idx=intersect(find(isfinite(ratio)),find(ratio)); 
nonzeroRatio = ratio(idx);
sum(nonzeroRatio<=1)
sum(nonzeroRatio>1)
title(sprintf('2015-11-09-3 %0.0f cells of %0.0f tested (%0.1f percent) are activated in the ''safe zone''',...
    sum(nonzeroRatio<=1),numel(nonzeroRatio),sum(nonzeroRatio<=1)/numel(nonzeroRatio)*100)); 

figure; 
plot(1:length(nonzeroRatio), cumsum(nonzeroRatio>1),'LineWidth',2,'Color','r')
hold on
plot(1:length(nonzeroRatio), (1:length(nonzeroRatio))/2, 'k')
plot(1:length(nonzeroRatio), cumsum(nonzeroRatio<1),'LineWidth',2,'Color',[0.5 0.9 0.5])
xlabel('hot spotty index')
ylabel('cumulative distribution of ratios');

figure; 
scatter(sortedHotSpotIdx(1:length(ratio)),ratio,'kd'); 
xlabel('Hot value sortedElectrodes(1)/sum(sortedElectrodes(2:8)) (higher -> more hot-spotty)');
ylabel('ratio: cell threshold : axon bundle threshold');
title(sprintf('2015-11-09-3 %0.0f cells of %0.0f tested (%0.1f percent) are activated in the ''safe zone''',...
    sum(nonzeroRatio<=1),numel(nonzeroRatio),sum(nonzeroRatio<=1)/numel(nonzeroRatio)*100)); 
hold on; line(get(gca,'XLim'),[1 1],'Color','r','LineStyle','--'); 
%%
% What are their activation thresholds at the maximum electrode?
figure; bar(allthresholds)
title('activation thresholds for top 10 hot-spotty cells'); 

playMovie512arrayAfterStimPattern_dots(['/Volumes/Analysis/2015-11-09-3/data001-data002/'], allmaxelecs(4),'movieindex',21) % No bundle activation, but there are delayed onset 'fireworks' that happen, indicating that a single axon is activated somewhere.
lastAmpBeforeBundle =[0.56 0.68 0.68 0.68 0.99 0.88 1.1 0.99 0.88 0.68]; 
figure; bar([allthresholds; lastAmpBeforeBundle]'); 
legend('activation threshold','bundle safe amplitude'); 
%% Look at a group of ON parasols from 2012-09-24-3
[actThresh,bundleDiff,cellToBundleRatio] = getThresholds201209243data008(); 
neurons = [2 33 168 197 303 332 436 601 647 6769 ...
    6936 6995 7008 7173 7218 7474 7507];

genActThreshSpatialMaps('/Volumes/Analysis/2012-09-24-3/data008/', neurons)

eipath = '/Volumes/Analysis/2012-09-24-3/data007/data007.ei';
[eiM,nIdList] = convertEiFileToMatrix(eipath);
eiAmps = squeeze(max(eiM,[],2) - min(eiM,[],2)); 
valsAllCells = zeros(size(eiM,1),size(eiM,3)); 
%%
[xc,yc] = getElectrodeCoords512(); 
customcolormap = summer(2);
for c = 1:size(neurons,2)
    nIdx = find(neurons(c) == nIdList); 
    ei = eiAmps(:,nIdx); 
    eiIdx =find(ei>2); 
    [mins,idx] = sort(min(eiM(:,:,nIdx),[],2)); 
    figure; 
    subplot(2,2,1); scatter(xc(eiIdx),yc(eiIdx),ei(eiIdx)*10,'filled'); axis image; axis off; 
    hold on; stimIdx = find(isfinite(actThresh(c,:))); 
    scatter(xc(stimIdx),yc(stimIdx),ei(stimIdx)*5,actThresh(c,stimIdx),'filled'); colorbar; caxis([0.5 3]);  
    [axon_x, axon_y, ~, soma_x, soma_y] = weighted_axon_poly_reg(ei);
    plot(axon_x, axon_y, 'Color', 'k','LineWidth',2);
    scatter(soma_x,soma_y,max(ei)+0.1,'k','filled');
    title(['n' num2str(neurons(c))]); 
    for i = 1:4
        text(xc(idx(i)),yc(idx(i)),num2str(i),'HorizontalAlignment','center'); 
    end
    subplot(2,2,3); scatter(xc(eiIdx),yc(eiIdx),ei(eiIdx)*10,'filled'); axis image; axis off; 
    hold on; bunIdx = find(isfinite(bundleDiff(c,:))); 
    scatter(xc(bunIdx),yc(bunIdx),ei(bunIdx)*5,bundleDiff(c,bunIdx),'filled'); colorbar; caxis([-3 3]);  colormap(customcolormap); 
    [axon_x, axon_y, ~, soma_x, soma_y] = weighted_axon_poly_reg(ei);
    plot(axon_x, axon_y, 'Color', 'k','LineWidth',2);
    scatter(soma_x,soma_y,max(ei)+0.1,'k','filled');
    title(['n' num2str(neurons(c)) ' bundle difference POSITIVE, GOOD']); 
   
    scaledVals = mins/mins(1); 
    subplot(2,2,2); bar(abs(mins(1:10))); xlabel('electrode order of signal strength'); 
    ylabel('measured voltage'); title('hot-spottiness (absolute signal)');
     subplot(2,2,4); bar(abs(scaledVals(1:10))); xlabel('electrode order of signal strength'); 
    ylabel('measured voltage'); title('hot-spottiness (percentage)');
    % Penalize edge cases.
    if abs(xc(idx(1)))>=915 || abs(yc(idx(1))) == max(yc)
        scaledVals(2) = 1; 
    end
    valsAllCells(:,c) = scaledVals; 
end


% Plot the ratio of cell thresholds to axon bundle thresholds for each cell
% at the electrode with maximum signal. 
for c = 1:size(neurons,2)
    nIdx = find(neurons(c) == nIdList); 
    ei = eiAmps(:,nIdx); 
    eiIdx =find(ei>2); 
    [mins,idx] = sort(min(eiM(:,:,nIdx),[],2)); 
    hotElectrode = idx(1); 
    hotSpotValue(c) = mins(1)/sum(mins(2:8)); 
    ratios(c) = cellToBundleRatio(c,hotElectrode);
end
figure; scatter(hotSpotValue,ratios,'kd'); 
xlabel('Hot value sortedElectrodes(1)/sum(sortedElectrodes(2:8)) (higher -> more hot-spotty)');
ylabel('ratio: cell threshold : axon bundle threshold');
title(sprintf('2012-09-24-3 %0.0f cells of %0.0f tested (%0.1f percent) are activated in the ''safe zone''',...
    sum(ratios<=1),numel(ratios),sum(ratios<=1)/numel(ratios)*100)); 
hold on; line(get(gca,'XLim'),[1 1],'Color','r','LineStyle','--'); 
%% Repeat for experiment 2015-10-06-3

fPath = '/Volumes/Analysis/2015-10-06-3/data001-data002-autosort-template000/';
eipath = '/Volumes/Analysis/2015-10-06-3/data000/data000.ei';
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/hand_sorted_bundle_thresholds/axonBundleThresholds_2015_10_06_3.mat')
bundleThresholds = axonBundleThresholds_2015_10_6_3; 
bundleThresholds(find(bundleThresholds==0)) = nan;
[eiM,nIdList] = convertEiFileToMatrix(eipath);
valsAllCells = zeros(size(eiM,1),size(eiM,3)); 

[xc,yc] = getElectrodeCoords512(); 
for c = 1:size(eiM,3)
    [mins,idx] = sort(min(eiM(:,:,c),[],2)); 
    scaledVals = mins/mins(1); 
    % Penalize edge cases.
    if abs(xc(idx(1)))>=915 || abs(yc(idx(1))) == max(yc)
        scaledVals(2) = 1; 
    end
    valsAllCells(:,c) = scaledVals; 
end

secondaryVals = sum(valsAllCells(2:3,:),1); 
[V,I] = sort(secondaryVals); 
sortedHotSpotIdx = 1./V; 
% Look at the first 10 neurons to see their 'hot spots'
nIdList(I(1:10))
eiAmps = squeeze(max(eiM,[],2) - min(eiM,[],2)); 
% figure; scatter(xc,yc,10,'k','filled'); axis image; axis off; 

% Calculate the thresholds of all the cells ranked by hot-spottiness
allthresholds= zeros(1,size(eiAmps,2)); 
allmaxelecs = zeros(1,size(eiAmps,2)); 
for ii = 1:size(eiAmps,2)
    idxx = find(eiAmps(:,I(ii))>1);
    [resort,idx]  = sort(min(eiM(:,:,I(ii)),[],2)); 
    fname = ['elecRespAuto_n' num2str( nIdList(I(ii))) '_p' num2str(idx(1)) '.mat']; 
    try
        temp = load(fullfile(fPath,fname));
        [threshold,~, ~]  = getThresholdFromAuto(temp.elecRespAuto);
        if threshold>5; threshold = 5; end
        if threshold<0; threshold = 5; end
        allthresholds(ii) = threshold;
        allmaxelecs(ii) = idx(1);
    catch
        disp([fname ' does not exist; signal of this neuron is ' num2str(resort(1))]);
        allthresholds(ii) = nan;
        allmaxelecs(ii) = nan;
    end   
end

figure; plot(allthresholds)
figure; hold on; 
colors = jet(size(allthresholds,2));
for ii = 1:size(allthresholds,2)
    if isfinite(allthresholds(ii))
        ratio(ii) = allthresholds(ii)/(bundleThresholds(allmaxelecs(ii))/1.1); 
        scatter(allthresholds(ii),bundleThresholds(allmaxelecs(ii))/1.1,50,sortedHotSpotIdx(ii),'filled')
    end
end
line(0:5,0:5,'Color','k');
xlabel('Somatic 50% activation threshold at electrode with max signal'); 
ylabel('Bundle safe zone,uA  (largest current with no bundle activation)'); 
c = colorbar;  colormap jet; 
c.Label.String = 'hot-spot value'; 
title('2015-10-06-3'); 
figure; hold on;
for ii = 1:size(allthresholds,2)
    if isfinite(allthresholds(ii))
        ratio(ii) = allthresholds(ii)/(bundleThresholds(allmaxelecs(ii))/1.1); 
        if isfinite(ratio(ii))
        if ratio(ii) > 1
            a=area([ii-1 ii],ratio(ii)*[1 1]); 
            a.FaceColor = [1 0 0]; 
             a.EdgeColor = [1 0 0];
        elseif ratio(ii)<1
             a=area([ii-1 ii],ratio(ii)*[1 1]); 
            a.FaceColor = [0.5 0.9 0.5]; 
            a.EdgeColor = [0.5 0.9 0.5]; 
        end
        end
    end
end
xlim([1 size(allthresholds,2)]); 
xlabel('cell hot-spottiness ranking, 1 is highest')
ylabel('ratio: cell threshold : axon bundle threshold')
 
line(get(gca,'XLim'),[1 1],'Color','k','LineStyle','--')
idx=intersect(find(isfinite(ratio)),find(ratio)); 
nonzeroRatio = ratio(idx);
sum(nonzeroRatio<=1)
sum(nonzeroRatio>1)
title(sprintf('2015-10-06-3 %0.0f cells of %0.0f tested (%0.1f percent) are activated in the ''safe zone''',...
    sum(nonzeroRatio<=1),numel(nonzeroRatio),sum(nonzeroRatio<=1)/numel(nonzeroRatio)*100)); 
figure; 
plot(1:length(nonzeroRatio), cumsum(nonzeroRatio>1),'LineWidth',2,'Color','r')
hold on
plot(1:length(nonzeroRatio), (1:length(nonzeroRatio))/2, 'k')
plot(1:length(nonzeroRatio), cumsum(nonzeroRatio<1),'LineWidth',2,'Color',[0.5 0.9 0.5])
xlabel('hot spotty index')
ylabel('cumulative distribution of ratios');

figure; 
plot(1:length(nonzeroRatio), cumsum(nonzeroRatio>1)./cumsum(nonzeroRatio<1),'LineWidth',2,'Color','r')
hold on
plot(1:length(nonzeroRatio), (1:length(nonzeroRatio))/2, 'k')
% plot(1:length(nonzeroRatio), cumsum(nonzeroRatio<1),'LineWidth',2,'Color',[0.5 0.9 0.5])
xlabel('hot spot rank')
ylabel('cumulative distribution of ratios');

figure; 
scatter(sortedHotSpotIdx(1:length(ratio)),ratio,'kd'); 
xlabel('Hot value sortedElectrodes(1)/sum(sortedElectrodes(2:8)) (higher -> more hot-spotty)');
ylabel('ratio: cell threshold : axon bundle threshold');
title(sprintf('2015-10-06-3 %0.0f cells of %0.0f tested (%0.1f percent) are activated in the ''safe zone''',...
    sum(nonzeroRatio<=1),numel(nonzeroRatio),sum(nonzeroRatio<=1)/numel(nonzeroRatio)*100)); 
hold on; line(get(gca,'XLim'),[1 1],'Color','r','LineStyle','--'); 