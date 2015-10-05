function[idx nlog1 P] = clustering_analysis_plots(X, km,gmm, nclusters, threed, twod, datarun000, cellids, tcnormnorm, mored, vec)

if(gmm)
    options = statset('Display','final'); %just fills in what options to be displayed by the output
    obj = gmdistribution.fit(X,nclusters,'Options',options, 'Replicates', 1, 'CovType', 'full', 'Start', vec);
    %obj = gmdistribution.fit(X,nclusters,'Options',options, 'Replicates', 10, 'CovType', 'full');
    [idx nlog1 P] = cluster(obj,X);
end
if(km)
    [idx, ctrs] = kmeans(X,nclusters, 'distance', 'sqEuclidean', 'display', 'final', 'replicates', 10, 'start', 'sample');
end

if(threed)
    figure
    for i = 1:nclusters
        if i ==1
            scatter3(X(idx==i,1),X(idx==i,2),X(idx==i,3),80,[0 0 1]);
        elseif i == 2;
            scatter3(X(idx==i,1),X(idx==i,2),X(idx==i,3),80,[1 0 0]);
        end
    hold on
    end
    legend('Cluster 1','Cluster 2', 'Cluster 3','Cluster 4','Cluster 5','Cluster 6', 'Cluster 7', 'Location','NW')
end

if(twod)
    figure
    for i = 1:nclusters
        if i ==1
            scatter(X(idx==i,1),X(idx==i,2),80,[0 0 1]);
        elseif i ==2
                scatter(X(idx==i,1),X(idx==i,2),80,[1 0 0]);
        end
    hold on
    end
    legend('Cluster 1','Cluster 2', 'Cluster 3','Cluster 4','Cluster 5','Location','NW')
end

if(mored)
end


%     figure(2)
%     for i = 1:nclusters
%     subplot(2,nclusters,i);
%     plot_rf_summaries(datarun000, cellids(idx==i), 'label', true)
%     subplot(2,nclusters,i+nclusters);
%     plot(tcnormnorm(:,idx==i), 'b');  
%     end
end



