function m = add(m, e, cluster)
%ADD Add an electrode to the model file.
%   Provide output from clustering algorithm
%   and all the acceptable clusters present on electrode e will
%   be written to the neurons file.
%
%   tamachado@salk.edu 4/3/08

if nargin == 1
    m.nObject.closeModel();
    return;
end

% We are only saving out information about the MAP, which is
% by design, always index 1 in each data structure.

mapIndex = cluster.mapIndex;
nClusters   = cluster.nClusters(1);

% Get means for clusters on electrode e
oMean = cluster.meanRecord(mapIndex);
oMean = oMean{1}';

% Get covariances for clusters on electrode e
oCov = cluster.covRecord(mapIndex);
oCov = oCov{1};

% Use means to figure out dimensionality of data
nDimensions = size(oMean,2);

% Weight each cluster equally if no weight information is provided (this is
% wrong and should be avoided)
if isfield(cluster, 'weights') == 0
    disp('Warning: Weight information not provided ... weighting all clusters equally')
    oProb = (1/nClusters) * ones(nClusters, 1);
else
    oProb = cluster.weights;
end

% Add up to 15 clusters. The largest clusters are added first.
[nil, sortedClusters] = sort(cluster.nSpikes, 'descend');

if nClusters > 15
    nClusters = 15;
end

prob = zeros(1,nClusters);
mean = zeros(nClusters, nDimensions);

for c = 1:nClusters
    % Only store the diagonal of each covariance matrix for each cluster
    temp = diag(oCov(:,:,sortedClusters(c)));
    
    % Reorganize mean and probability weight matrices
    prob(c) = oProb(sortedClusters(c));
    mean(c,:) = oMean(sortedClusters(c),:);
    
    for d = 1:nDimensions
        cov(d, c) = temp(d);
    end
end

cov = cov';

m.nObject.addElectrode(prob, mean, cov, e, nClusters);

disp(sprintf('%d clusters on electrode %d added to model file.\n', nClusters, e))


