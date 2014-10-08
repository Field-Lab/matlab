

%% parameters

pcWindow = [10 40]; %window of time in data to use for PCA (in samples)
datasetNo = 9;


%% 

d = getDatasetDetails(datasetNo);

cd(d.savePath)

mkdir dimensionality_reduction
cd dimensionality_reduction

nPatterns = length(d.patternNumbers);
nMovies = length(d.movieNumbers);

for i = 1:nPatterns
    for j = 1:nMovies
        
        %returns data traces: nTraces x nElectrodes x nSamples
        dataTraces=NS_ReadPreprocessedData(d.dataPath, '', 0, d.patternNumbers(i), d.movieNumbers(j));
        
        stimAmp = max(abs(getStimAmps(d.dataPath, d.patternNumbers(i), d.movieNumbers(i))));
        
        % center electrode only: finding PCs and explained variance
        data1Elec = squeeze(dataTraces(:, d.centerChannel, pcWindow(1):pcWindow(2)));
        size(data1Elec)
        data1Cov = cov(data1Elec);
        
        [coeff1, score1] = princomp(data1Elec);
        [coeff1_check, latent, explained1] = pcacov(data1Cov);
        
        runningTotal = 0;
        explained1Cumu = zeros(length(explained1));
        for k = 1:length(explained1)
            runningTotal = explained1(k) + runningTotal;
            explained1Cumu(k) = runningTotal;
        end
        
        % center and surrounding electrodes: finding PCs and explained variance
        electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
        channelsToUse=electrodeMap.getAdjacentsTo(d.centerChannel, 1)';
        
        temp = squeeze(dataTraces(:, channelsToUse, pcWindow(1):pcWindow(2)));
        dataSElec = reshape(permute(temp, [1 3 2]), size(temp, 1), size(temp, 2)*size(temp, 3));
        size(dataSElec)
        dataSCov = cov(dataSElec);
        
        [coeffS, scoreS] = princomp(dataSCov);
        [coeff2_check, latent, explainedS] = pcacov(dataSCov);
        
        runningTotal = 0;
        explainedSCumu = zeros(length(explainedS));
        for k = 1:length(explainedS)
            runningTotal = explainedS(k) + runningTotal;
            explainedSCumu(k) = runningTotal;
        end
        
        % center electrode only: dendrogram
        
        dist1Elec = pdist(data1Elec);
        
        link1ElecComp = linkage(dist1Elec, 'complete');
        link1ElecSing = linkage(dist1Elec, 'single');
        
        % determining largest single cluster after each new link
        nLink = size(link1ElecComp,1);
        
        clusterSizeComp = zeros(1, 2*nLink + 1);
        clusterSizeSing = zeros(1, 2*nLink + 1);
        largestClusterComp = zeros(1, nLink);
        largestClusterSing = zeros(1, nLink);
        dLargClusterComp = zeros(1, nLink - 1);
        dLargClusterSing = zeros(1, nLink - 1);
        clusterSizeComp(1:nLink+1) = ones(1, nLink+1);
        clusterSizeSing(1:nLink+1) = ones(1, nLink+1);
        linkSizeComp = zeros(1, nLink);
        linkSizeSing = zeros(1, nLink);
        
        for k = 1:nLink
            clusterSizeComp(100+k) = clusterSizeComp(link1ElecComp(k,1))+clusterSizeComp(link1ElecComp(k,2));
            clusterSizeSing(100+k) = clusterSizeSing(link1ElecSing(k,1))+clusterSizeSing(link1ElecSing(k,2));
            largestClusterComp(k) = max(clusterSizeComp);
            largestClusterSing(k) = max(clusterSizeSing);
            linkSizeComp(k) = link1ElecComp(k,3);
            linkSizeSing(k) = link1ElecSing(k,3);
            if k ~= 1
                dLargClusterComp(k-1) = (largestClusterComp(k)-largestClusterComp(k-1))/(linkSizeComp(k)-linkSizeComp(k-1));
                dLargClusterSing(k-1) = (largestClusterSing(k)-largestClusterSing(k-1))/(linkSizeSing(k)-linkSizeSing(k-1));
            end
        end
        
        
        figure
        subplot(2,1,1)
        hold on
        plot(linkSizeComp, largestClusterComp)
        plot(linkSizeComp(1:end-1), dLargClusterComp*(nLink/max(dLargClusterComp)),'r')
        hold off
        xlabel('link size')
        title('leaves in largest cluster (complete linkage)')
        
        subplot(2,1,2)
        hold on
        plot(linkSizeSing, largestClusterSing)
        plot(linkSizeSing(1:end-1), dLargClusterSing,'r')
        %plot(linkSizeSing(1:end-1), dLargClusterSing*(nLink/max(dLargClusterSing)),'r')
        hold off
        xlabel('link size')
        title('leaves in largest cluster (single linkage)')
        

        
        
        
%         linkClusters = cluster(link1Elec, 'cutoff', 0.7);
%         
%         dataLargestCluster = data1Elec(linkClusters == mode(linkClusters),:);
%         linkageLargestCluster = linkage(pdist(dataLargestCluster), 'complete');
%         c = max(linkageLargestCluster(:,3));

        nLinkBelowThreshComp = zeros(1, 500);
        nLinkBelowThreshSing = zeros(1, 500);
        threshValsComp = linspace(min(link1ElecComp(:,3)), max(link1ElecComp(:,3)), 500);
        threshValsSing = linspace(min(link1ElecSing(:,3)), max(link1ElecSing(:,3)), 500);
        
        for k = 1:500
            nLinkBelowThreshComp(k) = sum(squeeze(link1ElecComp(:,3)) <= threshValsComp(k));
            nLinkBelowThreshSing(k) = sum(squeeze(link1ElecSing(:,3)) <= threshValsSing(k));
        end
        
        figure
        subplot(2,1,1)
        hold on
        plot(threshValsComp, nLinkBelowThreshComp)
        plot(threshValsComp(2:end), diff(nLinkBelowThreshComp)*(nLink/max(diff(nLinkBelowThreshComp))))
        hold off
        title('links below threshold (complete linkage)')
        xlabel('threshold length')
        ylabel('number of links')
        
        subplot(2,1,2)
        hold on
        plot(threshValsSing, nLinkBelowThreshSing)
        plot(threshValsSing(2:end), diff(nLinkBelowThreshSing)*(nLink/max(diff(nLinkBelowThreshSing))), 'r')
        hold off
        title('links below threshold (single linkage)')
        xlabel('threshold length')
        ylabel('number of links')
        


        
        
%% center and surrounding electrodes: finding dendrogram
        
        distSElec = pdist(dataSElec);
        linkSElecComp = linkage(distSElec, 'complete');
        linkSElecSing = linkage(distSElec, 'single');
        
        
%% plotting
        
        % single electrode
        min1Elec = min(min(min(dataTraces(:, d.centerChannel, :))));
        max1Elec = max(max(max(dataTraces(:, d.centerChannel, :))));
        
        figure('Position', [100 100 800 1000])
%         uicontrol(gcf, 'Style', 'text', 'BackgroundColor', [1 1 1], 'String', [d.dataPath, 10,...
%             'pattern ', num2str(d.patternNumbers(i)), ', movie ',...
%             num2str(d.movieNumbers(j)), ', current amp: ', num2str(stimAmp)],...
%             'Position', [200 930 400 20])
        

        subplot(5,2,1)
        hold on
        plot(pcWindow(1):1:pcWindow(2), abs(coeff1(:,1))*(max1Elec-min1Elec)+max1Elec,'r-')
        plot(pcWindow(1):1:pcWindow(2), abs(coeff1(:,2))*(max1Elec-min1Elec)+max1Elec,'b-')
        for k = 1:size(dataTraces, 1);
            plot(squeeze(dataTraces(k, d.centerChannel, :)), 'k-')
        end
        % plots vertical lines between electrodes
        plot([pcWindow(1) pcWindow(1)], [min1Elec max1Elec], 'k--')
        plot([pcWindow(2) pcWindow(2)], [min1Elec max1Elec], 'k--')
        hold off
        title('raw data: primary electrode only')
        text(0, 1.5, [d.dataPath, 10,...
             'pattern ', num2str(d.patternNumbers(i)), ', movie ',...
             num2str(d.movieNumbers(j)), ', current amp: ', num2str(stimAmp)], 'Units','normalized')

        
        xlabel('samples')
        
        subplot(5,2,3)
        plot(score1(:,1), score1(:,2), 'k.')
        xlabel('PC1')
        ylabel('PC2')
        
        subplot(5,2,5)
        plot(explained1Cumu)
        set(gca, 'ylim', [0 100], 'xlim', [0 length(explained1)])
        ylabel('percent variance explained')
        xlabel('number of PCs')
        
        subplot(5,2,7)
        hDend1ElecComp = dendrogram(link1ElecComp, 0);
        title('complete linkage, primary electrode')
        set(gca, 'XTickMode', 'manual', 'XTickLabel',{''})
        ylabel('link distance')
        
        subplot(5,2,9)
        hDend1ElecSing = dendrogram(link1ElecSing, 0);
        title('single linkage, primary electrode')
        set(gca, 'XTickMode', 'manual', 'XTickLabel',{''})
        ylabel('link distance')
        
        %%%%%%%%%%%%%%%%% center + surrounding electrodes
        
        subplot(5,2,2)
        hold on
        plot(abs(coeffS(:,1))*(max1Elec-min1Elec)+max1Elec,'r-')
        plot(abs(coeffS(:,2))*(max1Elec-min1Elec)+max1Elec,'b-')
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
        xlabel('samples')
        
        
        subplot(5,2,4)
        plot(scoreS(:,1), scoreS(:,2), 'k.')
        xlabel('PC1')
        ylabel('PC2')
        
        
        subplot(5,2,6)
        plot(explainedSCumu)
        set(gca, 'ylim', [0 100], 'xlim', [0 length(explained1)])
        xlabel('number of PCs')
        
        subplot(5,2,8)
        hDendSElecComp = dendrogram(linkSElecComp, 0);
        title('complete linkage, primary and neighboring electrodes')
        set(gca, 'XTickMode', 'manual', 'XTickLabel',{''})
        ylabel('link distance')
        
        
        subplot(5,2,10)
        hDendSElecSing = dendrogram(linkSElecSing, 0);
        title('single linkage, primary and neighboring electrodes')
        set(gca, 'XTickMode', 'manual', 'XTickLabel',{''})
        ylabel('link distance')
        
        set(gcf,'PaperUnits','centimeters')
        xSize = 20; ySize=26;
        xLeft = (22-xSize)/2; yTop = (26-ySize)/2;
        set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
        
        saveas(gcf, ['p' num2str(d.patternNumbers(i)) 'm' num2str(d.movieNumbers(j)) '_dim_reduc.pdf'])
        close all
    end
end







