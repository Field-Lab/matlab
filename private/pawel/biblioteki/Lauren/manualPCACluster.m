function clusterIndex = manualPCACluster(traces)

% This function is designed to run PCA on the "traces" (after concatenating the traces
% on all electrodes for each pulse), and then hand-cluster the data.  After each cluster is
% selected, these traces are removed from the subset and PCA is performed again.  These iterations
% stop when all traces have been clustered.  The resulting array "clusterIndex" defines which traces
% belong to the first selected cluster (value of 1), which belong to the second selected cluster
% (value of 2), etc.
%
% arguments
%   traces: a 3D array, with traces X electrodes X samples
%
% returns
%   clusterIndex: a 1D array of integars representing which cluster corresponding traces belong to,
%   where clusters are numbered 1, 2, 3... nClusters
%
% written by Lauren Hruby, 2008-09-05
% last edited 2008-09-09

nTraces = size(traces, 1);
nElectrodes = size(traces, 2);
nSamples = size(traces, 3);

clusterIndex = zeros(1, nTraces);

% concatenates traces on different electrodes for each pulse
prinCompArray = zeros(nTraces, nElectrodes*nSamples);
for i = 1:nElectrodes
    prinCompArray(:, (i-1)*nSamples + 1 : i*nSamples) = traces(:,i,:);
end
s1=size(prinCompArray)
prinCompArray=traces;

clusterNumber = 1;
remainingTraces = nTraces;

while remainingTraces > 0
    
    [PCACoef, PCAScore] = princomp(prinCompArray);

    %displays 6 PC plots to choose from
    [PCx PCy] = PCChooser(PCAScore);
    
    if (PCx == 0 && PCy == 0) % if "cluster all remaining" was chosen
        index = 1:remainingTraces;
    else
        %displays data and chosen PC plot for hand-clustering
        color = cool(remainingTraces);
        figure; set(gcf,'position',[50,50,800,1200]); title('lasso a cluster'); subplot(2,1,1); hold on
        for j = 1:remainingTraces
            current = plot(prinCompArray(j,:));
            set(findobj(current,'Type','line'),'Color',color(j,:))
        end
        hold off; subplot(2,1,2)
        [selx, sely, index] = lasso(PCAScore(:,PCx), PCAScore(:,PCy));
        close; close
    end
    
    i = 1;
    j = 0;
    currentIndex = 1;
    while currentIndex <= length(index)
        if clusterIndex(i) == 0 %if trace was included in the current iteration of PCA/clustering
            j = j + 1;
            if index(currentIndex) == j %if the current trace is in the current cluster
                clusterIndex(i) = clusterNumber;
                prinCompArray(j, :) = prinCompArray(j, :)*0;
                currentIndex = currentIndex + 1;
            end
        end
        i = i + 1;
    end
    
    %removes the zeros rows from prinCompArray
    prinCompArray(~any(prinCompArray,2),:) = [];
    
    %updates cluster number and number of remaining traces
    clusterNumber = clusterNumber + 1;
    remainingTraces = sum(clusterIndex == 0);
end

end