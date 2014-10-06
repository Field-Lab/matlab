clear all
close all

%% parameters

pcWindow = [10 40]; %window of time in data to use for PCA (in samples)
datasetNo = 8;


%% 

d = getDatasetDetails(datasetNo);

cd(d.savePath)

mkdir linkage
cd linkage

nPatterns = length(d.patternNumbers);
nMovies = length(d.movieNumbers);

nDivisions = 100;
linkThresh1Elec = 20;
linkThreshSElec = 60;

for i = 1:nPatterns
    for j = 1:nMovies
        
        %returns data traces: nTraces x nElectrodes x nSamples
        %dataTraces=NS_ReadPreprocessedData(d.dataPath, '', 0, d.patternNumbers(i), d.movieNumbers(j));
        dataTraces=NS_ReadPreprocessedData(d.dataPath, '', 0, d.patternNumbers(i), 80);
        nTraces = size(dataTraces, 1);
        
        stimAmp = max(abs(getStimAmps(d.dataPath, d.patternNumbers(i), d.movieNumbers(j))));
        data1Elec = squeeze(dataTraces(:, d.centerChannel, pcWindow(1):pcWindow(2)));
        
        electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
        channelsToUse=electrodeMap.getAdjacentsTo(d.centerChannel, 1)';
        temp = squeeze(dataTraces(:, channelsToUse, pcWindow(1):pcWindow(2)));
        dataSElec = reshape(permute(temp, [1 3 2]), size(temp, 1), size(temp, 2)*size(temp, 3));
        
        % dendrogram
        
        dist1Elec = pdist(data1Elec);
        link1Elec = linkage(dist1Elec, 'single');
        
        distSElec = pdist(dataSElec);
        linkSElec = linkage(distSElec, 'single');
        
        % determining largest single cluster after each new link
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
        
        
        
        % finds leaves in 2 largest clusters at a particular link length threshold
        
        mainBranches1Elec = findSubThreshBranches(link1Elec, linkThresh1Elec);
        mainBranchesSElec = findSubThreshBranches(linkSElec, linkThreshSElec);
        
        
        % linear interpolation to make the x-values uniformly spaced
        
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
        
        
        % low-pass filter
        maFilter = fspecial('average',[3,1]);
        gFilter  = fspecial('gaussian',[3,1],1);
        
        largestClusterMA1Elec = conv(maFilter, largestClusterLin1Elec);
        largestClusterG1Elec  = conv(gFilter,  largestClusterLin1Elec);        
        largestClusterMASElec = conv(maFilter, largestClusterLinSElec);
        largestClusterGSElec  = conv(gFilter,  largestClusterLinSElec);
        
        largestClusterMA1Elec = largestClusterMA1Elec(2:end-1);
        largestClusterG1Elec  = largestClusterG1Elec(2:end-1);        
        largestClusterMASElec = largestClusterMASElec(2:end-1);
        largestClusterGSElec  = largestClusterGSElec(2:end-1);
        
        % differentiating
        
        dLargClustMA1Elec = diff(largestClusterMA1Elec);
        dLargClustG1Elec  = diff(largestClusterG1Elec);
        dLargClustMASElec = diff(largestClusterMASElec);
        dLargClustGSElec  = diff(largestClusterGSElec);
        dLargClustNF1Elec = diff(largestClusterLin1Elec);
        dLargClustNFSElec  = diff(largestClusterLinSElec);
        
        

        % determining number of links below thresholds

        nLinkBelowThresh1Elec = zeros(1, nDivisions);
        nLinkBelowThreshSElec = zeros(1, nDivisions);
        
        for k = 1:nDivisions
            nLinkBelowThresh1Elec(k) = sum(squeeze(link1Elec(:,3)) <= linkVals1Elec(k));
            nLinkBelowThreshSElec(k) = sum(squeeze(linkSElec(:,3)) <= linkValsSElec(k));
        end
        

        
 
        
        
%% plotting
        
        % single electrode
        min1Elec = min(min(min(dataTraces(:, d.centerChannel, :))));
        max1Elec = max(max(max(dataTraces(:, d.centerChannel, :))));
        
        %figure('Position', [100 100 800 1000], 'Visible', 'off') 
        figure('Position', [100 100 800 1000])
        subplot(5,2,1)
        
        text(-0.2, 1.8, [d.dataPath, 10,...
            'pattern ', num2str(d.patternNumbers(i)), ', movie ',...
            num2str(d.movieNumbers(j)), ', current amp: ', num2str(stimAmp)], 'Units','normalized')
        
        
        hold on
        for k = 1:size(dataTraces, 1);
            plot(squeeze(dataTraces(k, d.centerChannel, :)), 'k-')
        end
        % plots vertical lines around analyzed region
        plot([pcWindow(1) pcWindow(1)], [min1Elec max1Elec], 'k--')
        plot([pcWindow(2) pcWindow(2)], [min1Elec max1Elec], 'k--')
        hold off
        title('raw data: primary electrode only')
        xlabel('samples', 'FontSize', 8)  

        
        subplot(5,2,3)
        hDend1ElecSing = dendrogram(link1Elec, 0, 'colorthreshold', linkThresh1Elec);
        title('single linkage, primary electrode')
        set(gca, 'XTickMode', 'manual', 'XTickLabel',{''})
        ylabel('link distance', 'FontSize', 8)
        
        
        subplot(5,2,5)
        hold on
        plot(linkVals1Elec, nLinkBelowThresh1Elec, 'k')
        plot(linkVals1Elec(2:end), diff(nLinkBelowThresh1Elec)*(nLink/max(diff(nLinkBelowThresh1Elec))))
        hold off
        title('links below threshold')
        xlabel('threshold length', 'FontSize', 8)
        ylabel('number of links', 'FontSize', 8)
        
        subplot(5,2,7)
        hold on
        plot(linkSize1Elec, largestCluster1Elec, 'k')
        plot(linkVals1Elec(2:end-2), dLargClustG1Elec(1:end-2)*max(largestCluster1Elec)/max(dLargClustNF1Elec),'b')
        plot(linkVals1Elec(2:end-2), dLargClustNF1Elec(1:end-2)*max(largestCluster1Elec)/max(dLargClustNF1Elec),'m')
        hold off
        xlabel('link size', 'FontSize', 8)
        ylabel('cluster size or scaled difference', 'FontSize', 8)
        title('leaves in largest cluster')
        

        subplot(5,2,9)
        hold on
        for k = 1:size(dataTraces, 1);
            %if sum(leafIDInCluster1Elec{mainClusterID1Elec} == k)
            if sum(mainBranches1Elec{1} == k)
                plot(squeeze(dataTraces(k, d.centerChannel, :)), 'k-')
            end
        end
        for k = 1:size(dataTraces, 1);
            %if sum(leafIDInCluster1Elec{secondClusterID1Elec} == k)
            if sum(mainBranches1Elec{2} == k)
                plot(squeeze(dataTraces(k, d.centerChannel, :)), 'm-')
            end
        end
        % plots vertical lines around analyzed region
        plot([pcWindow(1) pcWindow(1)], [min1Elec max1Elec], 'k--')
        plot([pcWindow(2) pcWindow(2)], [min1Elec max1Elec], 'k--')
        hold off
        title(['traces from two largest clusters', 10, 'below threshold link size = ', num2str(linkThresh1Elec)])
        
        
        
        %%%%%%%%%%%%%%%%% center + surrounding electrodes
        
        subplot(5,2,2)
        hold on
        for k = 1:size(dataTraces, 1);
            plot(squeeze(dataSElec(k, :)), 'k-')
        end
        
        % plots vertical lines between electrodes
        for k = 1:length(channelsToUse-1)
            position = (pcWindow(2)-pcWindow(1)+1)*k;
            plot([position, position], [min1Elec max1Elec], 'k--')
            plot([position+1, position+1], [min1Elec max1Elec], 'k--')
        end
        
        hold off
        title('raw data: primary and neighboring electrodes')
        xlabel('samples', 'FontSize', 8)

        
        subplot(5,2,4)
        hDendSElec = dendrogram(linkSElec, 0, 'colorthreshold', linkThreshSElec);
        title('single linkage, primary and neighboring electrodes')
        set(gca, 'XTickMode', 'manual', 'XTickLabel',{''})
        ylabel('link distance', 'FontSize', 8)
     
        
        subplot(5,2,6)
        hold on
        plot(linkValsSElec, nLinkBelowThreshSElec, 'k')
        plot(linkValsSElec(2:end), diff(nLinkBelowThreshSElec)*(nLink/max(diff(nLinkBelowThreshSElec))))
        hold off
        title('links below threshold')
        xlabel('threshold length', 'FontSize', 8)
        ylabel('number of links', 'FontSize', 8)
        
        
        subplot(5,2,8)
        hold on
        plot(linkSizeSElec, largestClusterSElec, 'k')
        plot(linkValsSElec(2:end-2), dLargClustGSElec(1:end-2)*max(largestClusterSElec)/max(dLargClustNFSElec),'b')
        plot(linkValsSElec(2:end-2), dLargClustNFSElec(1:end-2)*max(largestClusterSElec)/max(dLargClustNFSElec),'m')
        hold off
        xlabel('link size', 'FontSize', 8)
        ylabel('cluster size or scaled difference', 'FontSize', 8)
        title('leaves in largest cluster')
        
        subplot(5,2,10)
        hold on
        for k = 1:size(dataTraces, 1);
            %if sum(leafIDInClusterSElec{mainClusterIDSElec} == k)
            if sum(mainBranchesSElec{1} == k)
                plot(squeeze(dataTraces(k, d.centerChannel, :)), 'k-')
            end
        end
        
        for k = 1:size(dataTraces, 1);
            %if sum(leafIDInClusterSElec{secondClusterIDSElec} == k)
            if sum(mainBranchesSElec{2} == k)
                plot(squeeze(dataTraces(k, d.centerChannel, :)), 'm-')
            end
        end
        % plots vertical lines around analyzed region
        plot([pcWindow(1) pcWindow(1)], [min1Elec max1Elec], 'k--')
        plot([pcWindow(2) pcWindow(2)], [min1Elec max1Elec], 'k--')
        hold off
        title(['traces from two largest clusters', 10, 'below threshold link size = ', num2str(linkThreshSElec)])
        
        keyboard
        
        set(gcf,'PaperUnits','centimeters')
        xSize = 20; ySize=26;
        xLeft = (22-xSize)/2; yTop = (26-ySize)/2;
        set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
        
        saveas(gcf, ['p' num2str(d.patternNumbers(i)) 'm' num2str(d.movieNumbers(j))...
            '_linkage_' num2str(linkThresh1Elec) '-' num2str(linkThreshSElec) '.pdf'])
        
        % keyboard
        close all
    end
end







