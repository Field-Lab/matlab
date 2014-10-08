clear all
close all

%% parameters

pcWindow = [10 40]; %window of time in data to use for PCA (in samples)


%% getting artifacts

dataPath1 = '/Volumes/Lee/Analysis/Lauren/';
dataPath2a = '2008-08-26-0';
dataPath2b = 'data006';
dataPath = [dataPath1 dataPath2a '/' dataPath2b];
artifactPath = '/Volumes/Lee/Analysis/Lauren/2008-08-26-0/data011';

savePath = '/netapp/snle/home/lhruby/Documents/projects/data_analysis/automated/analysis/artifactLinkages';

%clusterPath = '/Volumes/Lee/Analysis/Lauren/2008-11-10-3/clusters011artifacts';
clusterPath = '/Volumes/Lee/Analysis/Lauren/2008-08-26-0/clusters006artifacts';


centerChannel = 56;
patternNumber = 56;
PW = 200;

patternNumbers = patternNumber;
movieNumbers = 3:3:75;

artPatternNumbers = patternNumber;
%artMovieNumbers = 2:3:74;
artMovieNumbers = 3:3:75;



%%

cd(savePath)
mkdir([dataPath2a '-' dataPath2b])
cd([dataPath2a '-' dataPath2b])

nPatterns = length(patternNumbers);
nMovies = length(movieNumbers);

stimAmp = zeros(nPatterns, nMovies);

linkSize80Thresh = zeros(nPatterns, nMovies);
linkSize80Thresh1ElecArt = zeros(nPatterns, nMovies);
linkSize80ThreshSElecArt = zeros(nPatterns, nMovies);

nDivisions = 100;
nDivisionsArt = 50;
%linkThresh1Elec = 20;
%linkThreshSElec = 80;

for i = 1:nPatterns
    for j = 1:nMovies
        
        %returns data traces: nTraces x nElectrodes x nSamples
        dataTraces=NS_ReadPreprocessedData(dataPath, '', 0, patternNumbers(i), movieNumbers(j));
        artifactTraces = NS_ReadPreprocessedData(artifactPath, '', 0, artPatternNumbers(i), artMovieNumbers(j));
        clusterIndex = NS_ReadClusterFile(clusterPath, movieNumbers(j), patternNumbers(i));
        
        % removes dataTraces that have been clustered out as having spike waveforms in them
        dataTraces = dataTraces(clusterIndex==1, :, :);
        
        nTraces = size(dataTraces, 1);
        nArtifacts = size(artifactTraces, 1);
        
        stimAmp(i,j) = max(abs(getStimAmps(dataPath, patternNumbers(i), movieNumbers(j))));
        data1Elec = squeeze(dataTraces(:, centerChannel, pcWindow(1):pcWindow(2)));
        artifact1Elec = squeeze(artifactTraces(:, centerChannel, pcWindow(1):pcWindow(2)));
        
        electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
        channelsToUse=electrodeMap.getAdjacentsTo(centerChannel, 1)';
        temp = squeeze(dataTraces(:, channelsToUse, pcWindow(1):pcWindow(2)));
        dataSElec = reshape(permute(temp, [1 3 2]), size(temp, 1), size(temp, 2)*size(temp, 3));
        temp = squeeze(artifactTraces(:, channelsToUse, pcWindow(1):pcWindow(2)));
        artifactSElec = reshape(permute(temp, [1 3 2]), size(temp, 1), size(temp, 2)*size(temp, 3));

        
        % dendrogram
        
        dist1Elec = pdist(data1Elec);
        link1Elec = linkage(dist1Elec, 'single');
        dist1ElecArt = pdist(artifact1Elec);
        link1ElecArt = linkage(dist1ElecArt, 'single');
        
        distSElec = pdist(dataSElec);
        linkSElec = linkage(distSElec, 'single');
        distSElecArt = pdist(artifactSElec);
        linkSElecArt = linkage(distSElecArt, 'single');
        
        
        % determining largest single cluster after each new link (dataTraces)
        nLink = size(link1Elec,1);
       
        clusterSize1Elec = zeros(1, 2*nLink + 1);
        clusterSizeSElec = zeros(1, 2*nLink + 1);
        largestCluster1Elec = zeros(1, nLink);
        largestClusterSElec = zeros(1, nLink);
        clusterSize1Elec(1:nLink+1) = ones(1, nLink+1);
        clusterSizeSElec(1:nLink+1) = ones(1, nLink+1);
        linkSize1Elec = zeros(1, nLink);
        linkSizeSElec = zeros(1, nLink);
        leafIDInCluster1Elec = cell(1, 2*nLink+1);
        leafIDInClusterSElec = cell(1, 2*nLink+1);
        largestClusterID1Elec = zeros(1, nLink);
        largestClusterIDSElec = zeros(1, nLink);
        for k = 1:nLink
            leafIDInCluster1Elec{k} = k;
            leafIDInClusterSElec{k} = k;
        end
        
        for k = 1:nLink
            clusterSize1Elec(nTraces+k) = clusterSize1Elec(link1Elec(k,1))+clusterSize1Elec(link1Elec(k,2));
            clusterSizeSElec(nTraces+k) = clusterSizeSElec(linkSElec(k,1))+clusterSizeSElec(linkSElec(k,2));
            largestCluster1Elec(k) = max(clusterSize1Elec);
            largestClusterSElec(k) = max(clusterSizeSElec);
            largestClusterID1Elec(k) = find(clusterSize1Elec==max(clusterSize1Elec),1);
            largestClusterIDSElec(k) = find(clusterSizeSElec==max(clusterSizeSElec),1);
            linkSize1Elec(k) = link1Elec(k,3);
            linkSizeSElec(k) = linkSElec(k,3);
            leafIDInCluster1Elec{nTraces+k} = [leafIDInCluster1Elec{link1Elec(k,1)} leafIDInCluster1Elec{link1Elec(k,2)}];
            leafIDInClusterSElec{nTraces+k} = [leafIDInClusterSElec{linkSElec(k,1)} leafIDInClusterSElec{linkSElec(k,2)}];
        end
        
        
        % determining largest single cluster after each new link (TTX artifacts)
        nLinkArt = size(link1ElecArt,1);
       
        clusterSize1ElecArt = zeros(1, 2*nLinkArt + 1);
        clusterSizeSElecArt = zeros(1, 2*nLinkArt + 1);
        largestCluster1ElecArt = zeros(1, nLinkArt);
        largestClusterSElecArt = zeros(1, nLinkArt);
        clusterSize1ElecArt(1:nLinkArt+1) = ones(1, nLinkArt+1);
        clusterSizeSElecArt(1:nLinkArt+1) = ones(1, nLinkArt+1);
        linkSize1ElecArt = zeros(1, nLinkArt);
        linkSizeSElecArt = zeros(1, nLinkArt);
        leafIDInCluster1ElecArt = cell(1, 2*nLinkArt+1);
        leafIDInClusterSElecArt = cell(1, 2*nLinkArt+1);
        largestClusterID1ElecArt = zeros(1, nLinkArt);
        largestClusterIDSElecArt = zeros(1, nLinkArt);
        for k = 1:nLinkArt
            leafIDInCluster1ElecArt{k} = k;
            leafIDInClusterSElecArt{k} = k;
        end
        
        for k = 1:nLinkArt
            clusterSize1ElecArt(nArtifacts+k) = clusterSize1ElecArt(link1ElecArt(k,1))+clusterSize1ElecArt(link1ElecArt(k,2));
            clusterSizeSElecArt(nArtifacts+k) = clusterSizeSElecArt(linkSElecArt(k,1))+clusterSizeSElecArt(linkSElecArt(k,2));
            largestCluster1ElecArt(k) = max(clusterSize1ElecArt);
            largestClusterSElecArt(k) = max(clusterSizeSElecArt);
            largestClusterID1ElecArt(k) = find(clusterSize1ElecArt==max(clusterSize1ElecArt),1);
            largestClusterIDSElecArt(k) = find(clusterSizeSElecArt==max(clusterSizeSElecArt),1);
            linkSize1ElecArt(k) = link1ElecArt(k,3);
            linkSizeSElecArt(k) = linkSElecArt(k,3);
            leafIDInCluster1ElecArt{nArtifacts+k} = [leafIDInCluster1ElecArt{link1ElecArt(k,1)} leafIDInCluster1ElecArt{link1ElecArt(k,2)}];
            leafIDInClusterSElecArt{nArtifacts+k} = [leafIDInClusterSElecArt{linkSElecArt(k,1)} leafIDInClusterSElecArt{linkSElecArt(k,2)}];
        end
        

   
        % linear interpolation to make the x-values uniformly spaced (dataTraces)
        
        % removes any x-value repeats
        for k = nLink - 1 : -1 : 1
            if linkSize1Elec(k) == linkSize1Elec(k+1)
                linkSize1Elec(k+1) = [];
                largestCluster1Elec(k) = [];
            end
            if linkSizeSElec(k) == linkSizeSElec(k+1)
                linkSizeSElec(k+1) = [];
                largestClusterSElec(k) = [];
            end
        end
        
        linkVals1Elec = linspace(min(link1Elec(:,3)), max(link1Elec(:,3)), nDivisions);
        linkValsSElec = linspace(min(linkSElec(:,3)), max(linkSElec(:,3)), nDivisions);
        largestClusterLin1Elec = interp1(linkSize1Elec, largestCluster1Elec, linkVals1Elec);
        largestClusterLinSElec = interp1(linkSizeSElec, largestClusterSElec, linkValsSElec);
        
        
        % linear interpolation to make the x-values uniformly spaced (ttx artifacts)
        
        % removes any x-value repeats
        for k = nLinkArt - 1 : -1 : 1
            if linkSize1ElecArt(k) == linkSize1ElecArt(k+1)
                linkSize1ElecArt(k+1) = [];
                largestCluster1ElecArt(k) = [];
            end
            if linkSizeSElecArt(k) == linkSizeSElecArt(k+1)
                linkSizeSElecArt(k+1) = [];
                largestClusterSElecArt(k) = [];
            end
        end
        
        linkVals1ElecArt = linspace(min(link1ElecArt(:,3)), max(link1ElecArt(:,3)), nDivisionsArt);
        linkValsSElecArt = linspace(min(linkSElecArt(:,3)), max(linkSElecArt(:,3)), nDivisionsArt);
        largestClusterLin1ElecArt = interp1(linkSize1ElecArt, largestCluster1ElecArt, linkVals1ElecArt);
        largestClusterLinSElecArt = interp1(linkSizeSElecArt, largestClusterSElecArt, linkValsSElecArt);
        
        
        % low-pass filter
        gFilter = fspecial('gaussian',[5,1],1);
        
        largestClusterG1Elec  = conv(gFilter,  largestClusterLin1Elec);        
        largestClusterGSElec  = conv(gFilter,  largestClusterLinSElec);
              
        largestClusterG1ElecArt  = conv(gFilter,  largestClusterLin1ElecArt);        
        largestClusterGSElecArt  = conv(gFilter,  largestClusterLinSElecArt);
        
        largestClusterG1Elec  = largestClusterG1Elec(3:end-2);        
        largestClusterGSElec  = largestClusterGSElec(3:end-2);
        
        largestClusterG1ElecArt  = largestClusterG1ElecArt(3:end-2);        
        largestClusterGSElecArt  = largestClusterGSElecArt(3:end-2);
        
        % differentiating
        

        dLargClustG1Elec  = diff(largestClusterG1Elec);
        dLargClustGSElec  = diff(largestClusterGSElec);
        dLargClustNF1Elec = diff(largestClusterLin1Elec);
        dLargClustNFSElec = diff(largestClusterLinSElec);
        
        dLargClustG1ElecArt  = diff(largestClusterG1ElecArt);
        dLargClustGSElecArt  = diff(largestClusterGSElecArt);
        dLargClustNF1ElecArt = diff(largestClusterLin1ElecArt);
        dLargClustNFSElecArt = diff(largestClusterLinSElecArt);
        
       
        % determining number of links below thresholds

        nLinkBelowThresh1Elec = zeros(1, nDivisions);
        nLinkBelowThreshSElec = zeros(1, nDivisions);
        nLinkBelowThresh1ElecArt = zeros(1, nDivisionsArt);
        nLinkBelowThreshSElecArt = zeros(1, nDivisionsArt);
        
        threshFound = 0;
        for k = 1:nDivisions
            nLinkBelowThresh1Elec(k) = sum(squeeze(link1Elec(:,3)) <= linkVals1Elec(k));
            nLinkBelowThreshSElec(k) = sum(squeeze(linkSElec(:,3)) <= linkValsSElec(k));
            if nLinkBelowThresh1Elec(k) >= 0.8*nLink && ~threshFound
                linkSize80Thresh(i,j) = linkVals1Elec(k);
                threshFound = 1;
            end
        end
        
        threshFound1Elec = 0;
        threshFoundSElec = 0;
        for k = 1:nDivisionsArt
            nLinkBelowThresh1ElecArt(k) = sum(squeeze(link1ElecArt(:,3)) <= linkVals1ElecArt(k));
            nLinkBelowThreshSElecArt(k) = sum(squeeze(linkSElecArt(:,3)) <= linkValsSElecArt(k));
            if nLinkBelowThresh1ElecArt(k) >= 0.8*nLinkArt && ~threshFound1Elec
                linkSize80Thresh1ElecArt(i,j) = linkVals1ElecArt(k);
                threshFound1Elec = 1;
            end
            if nLinkBelowThreshSElecArt(k) >= 0.8*nLinkArt && ~threshFoundSElec
                linkSize80ThreshSElecArt(i,j) = linkValsSElecArt(k);
                threshFoundSElec = 1;
            end
        end
        
        
        
        % low-pass filter
        
        nLinkBelowThresh1ElecFilt = conv(gFilter, nLinkBelowThresh1Elec);
        nLinkBelowThreshSElecFilt = conv(gFilter, nLinkBelowThreshSElec);
        nLinkBelowThresh1ElecFiltArt = conv(gFilter, nLinkBelowThresh1ElecArt);
        nLinkBelowThreshSElecFiltArt = conv(gFilter, nLinkBelowThreshSElecArt);

        nLinkBelowThresh1ElecFilt = nLinkBelowThresh1ElecFilt(3:end-2);
        nLinkBelowThreshSElecFilt = nLinkBelowThreshSElecFilt(3:end-2);
        nLinkBelowThresh1ElecFiltArt = nLinkBelowThresh1ElecFiltArt(3:end-2);
        nLinkBelowThreshSElecFiltArt = nLinkBelowThreshSElecFiltArt(3:end-2);
        
        dLinkBelowThresh1Elec = diff(nLinkBelowThresh1ElecFilt);
        dLinkBelowThreshSElec = diff(nLinkBelowThreshSElecFilt);
        dLinkBelowThresh1ElecArt = diff(nLinkBelowThresh1ElecFiltArt);
        dLinkBelowThreshSElecArt = diff(nLinkBelowThreshSElecFiltArt);
        
        
%% plotting
        
        % single electrode
        min1Elec = min(min(min(dataTraces(:, centerChannel, :))));
        max1Elec = max(max(max(dataTraces(:, centerChannel, :))));
        
        figure('Position', [100 100 800 1000], 'visible', 'off') 

        subplot(5,2,1)
        
        text(-0.2, 1.5, [dataPath2a, '/', dataPath2b, 10,...
            'pattern ', num2str(patternNumbers(i)), ', movie ',...
            num2str(movieNumbers(j)), ', current amp: ', num2str(stimAmp(i,j)), ', PW: ', num2str(PW)], 'Units','normalized')
        
        
        hold on
        for k = 1:size(dataTraces, 1);
            hdata = plot(squeeze(dataTraces(k, centerChannel, :)), 'k-');
        end
        for k = 1:size(artifactTraces, 1);
            hart = plot(squeeze(artifactTraces(k, centerChannel, :)), 'm-');
        end
        % plots vertical lines around analyzed region
        plot([pcWindow(1) pcWindow(1)], [min1Elec max1Elec], 'k--')
        plot([pcWindow(2) pcWindow(2)], [min1Elec max1Elec], 'k--')
        hold off
        title('raw artifacts: primary electrode only')
        xlabel('samples', 'FontSize', 8)
        legendlabels{1} = 'no TTX'; %#ok<AGROW>
        legendlabels{2} = 'TTX'; %#ok<AGROW>
        legend([hdata hart], legendlabels, 'Location', 'NorthEast', 'FontSize', 8)

        
        subplot(5,2,3)
        h = dendrogram(link1Elec, 0);
        set(findobj(h,'Type','line'),'Color', [0 0 0])
        title('primary electrode, artifact-only data traces')
        set(gca, 'XTickMode', 'manual', 'XTickLabel',{''})
        ylabel('link distance', 'FontSize', 8)
        
        subplot(5,2,5)
        h = dendrogram(link1ElecArt, 0);
        set(findobj(h,'Type','line'),'Color', [0 0 0])
        title('primary electrode, ttx artifacts')
        set(gca, 'XTickMode', 'manual', 'XTickLabel',{''})
        ylabel('link distance', 'FontSize', 8)
        
        
        subplot(5,2,7)
        hold on
        plot(linkVals1Elec, nLinkBelowThresh1Elec*100/nLink, 'k', 'LineWidth', 2)
        plot(linkVals1Elec(2:end), dLinkBelowThresh1Elec*(100/max(dLinkBelowThresh1Elec)), 'k')
        plot(linkVals1ElecArt, nLinkBelowThresh1ElecArt*100/nLinkArt, 'm', 'LineWidth', 2)
        plot(linkVals1ElecArt(2:end), dLinkBelowThresh1ElecArt*(100/max(dLinkBelowThresh1ElecArt)), 'm')
        hold off
        axis tight
        set(gca, 'YLim', [0 110])
        title('links below threshold')
        xlabel('threshold length', 'FontSize', 8)
        ylabel('percent of links', 'FontSize', 8)
        
        subplot(5,2,9)
        hold on
        plot(linkSize1Elec, largestCluster1Elec*100/nLink, 'k', 'LineWidth', 2)
        plot(linkVals1Elec(2:end-2), dLargClustG1Elec(1:end-2)*100/max(dLargClustG1Elec),'k')
        plot(linkSize1ElecArt, largestCluster1ElecArt*100/nLinkArt, 'm', 'LineWidth', 2)
        plot(linkVals1ElecArt(2:end-2), dLargClustG1ElecArt(1:end-2)*100/max(dLargClustG1ElecArt),'m')
        hold off
        axis tight
        set(gca, 'YLim', [0 110])
        xlabel('link size', 'FontSize', 8)
        ylabel('cluster size, in % total leaves', 'FontSize', 8)
        title('leaves in largest cluster')
        
        
        
        
        %%%%%%%%%%%%%%%%% center + surrounding electrodes
        
        subplot(5,2,2)
        hold on
        for k = 1:size(dataTraces, 1);
            plot(squeeze(dataSElec(k, :)), 'k-')
        end
        for k = 1:size(artifactTraces, 1);
            plot(squeeze(artifactSElec(k, :)), 'm-')
        end
        
        % plots vertical lines between electrodes
        for k = 1:length(channelsToUse-1)
            position = (pcWindow(2)-pcWindow(1)+1)*k;
            plot([position, position], [min1Elec max1Elec], 'k--')
            plot([position+1, position+1], [min1Elec max1Elec], 'k--')
        end
        
        hold off
        title('raw artifacts: primary and neighboring electrodes')
        xlabel('samples', 'FontSize', 8)

        
        subplot(5,2,4)
        h = dendrogram(linkSElec, 0);
        set(findobj(h,'Type','line'),'Color', [0 0 0])
        title('primary and neighboring electrodes, artifact-only data traces')
        set(gca, 'XTickMode', 'manual', 'XTickLabel',{''})
        ylabel('link distance', 'FontSize', 8)
        
        subplot(5,2,6)
        h = dendrogram(linkSElecArt, 0);
        set(findobj(h,'Type','line'),'Color', [0 0 0])
        title('primary and neighboring electrodes, ttx artifacts')
        set(gca, 'XTickMode', 'manual', 'XTickLabel',{''})
        ylabel('link distance', 'FontSize', 8)
        
        
        subplot(5,2,8)
        hold on
        %plot(linkValsSElec, nLinkBelowThreshSElec*100/nLink, 'k')
        %plot(linkValsSElec(2:end), dLinkBelowThreshSElec*(100/max(dLinkBelowThreshSElec)), 'k')
        plot(linkValsSElecArt, nLinkBelowThreshSElecArt*100/nLinkArt, 'm', 'LineWidth', 2)
        plot(linkValsSElecArt(2:end), dLinkBelowThreshSElecArt*(100/max(dLinkBelowThreshSElecArt)), 'm')
        hold off
        axis tight
        set(gca, 'YLim', [0 110])
        title('links below threshold')
        xlabel('threshold length', 'FontSize', 8)
        ylabel('percent of links', 'FontSize', 8)
        
        subplot(5,2,10)
        hold on
        %plot(linkSizeSElec, largestClusterSElec*100/nLink, 'k')
        %plot(linkValsSElec(2:end-2), dLargClustGSElec(1:end-2)*100/max(dLargClustNFSElec),'k')
        plot(linkSizeSElecArt, largestClusterSElecArt*100/nLinkArt, 'm', 'LineWidth', 2)
        plot(linkValsSElecArt(2:end-2), dLargClustGSElecArt(1:end-2)*100/max(dLargClustGSElecArt),'m')
        hold off
        axis tight
        set(gca, 'YLim', [0 110])
        xlabel('link size', 'FontSize', 8)
        ylabel('cluster size, in % total leaves', 'FontSize', 8)
        title('leaves in largest cluster')
        
       
        
        
        set(gcf,'PaperUnits','centimeters')
        xSize = 20; ySize=26;
        xLeft = (22-xSize)/2; yTop = (26-ySize)/2;
        set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
        
        saveas(gcf, ['p' num2str(patternNumbers(i)) 'm' num2str(movieNumbers(j)) '_art_linkage.pdf'])
        
        %keyboard
        close all
    end
end

figure
hold on
plot(stimAmp, linkSize80Thresh, 'k')
plot(stimAmp, linkSize80Thresh1ElecArt, 'm')
plot(stimAmp, linkSize80ThreshSElecArt, 'c')
hold off

legend('single electrode, artifact-only data', 'single electrode, TTX artifacts', 'single + surrounding electrodes, TTX artifacts')
xlabel('stimulation current amplitude')
ylabel('link length threshold that includes 80+% of links')
title(['summary of link length threshold required to capture 80% of artifact links', 10,...
    dataPath2a, '/', dataPath2b, ' pattern ', num2str(patternNumbers(i)), ', PW:' num2str(PW), '\mus'])
set(gca, 'YLim', [0 100])

saveName = ['p' num2str(patternNumbers(i)) 'PW' num2str(PW)];
saveas(gcf, [saveName 'link_thresh_summary.pdf'])
dataSaveName = [dataPath2a, '-', dataPath2b, saveName, '_linkThresh80.mat'];
save(dataSaveName, 'stimAmp', 'linkSize80Thresh', 'linkSize80Thresh1ElecArt', 'linkSize80ThreshSElecArt')


