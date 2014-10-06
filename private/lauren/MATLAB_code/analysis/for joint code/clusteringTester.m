%% script to test clustering types/parameters

clear all
close all

DataPath = '/Analysis.noindex/Lauren/2008-08-26-0/data006_proba';
nPatterns = 5;
nMovies = 5;

dataTraces = cell(nPatterns, nMovies);
channels = cell(nPatterns, nMovies);
%clusterIndex = cell(nPatterns, nMovies);
%PCAScore = cell(nPatterns, nMovies);

%% clustering parameters


verbose = true; %specifies whether test output is displayed during clustering
iterations =  250; %Iteration settings (values from config.m)
DIMS = 5; %number of PCs included in clustering
nBurnIn = 125;

% CRP/IMM Hyperparameters (values from config.m)
% Use sample_igmm_prior.m to determine the best priors to use

% Parameters for the gamma function (sampled from when computing 
% alpha in the chinese restaurant process sampler)
a_0 = 1; %increasing this value increases the number of clusters predicted by priors 
b_0 = 1; %increasing this value increases the number of clusters predicted by priors

% Priors on cluster means
mu_0 = zeros(DIMS,1);
k_0 = 0.1;

% Priors on cluster covariances
lambda_0 = eye(DIMS)*1000;
v_0 = 10;




%% loading data
for i = 1:nPatterns
    if i == 9 || i == 25 || i == 57
    else
        disp(sprintf('loading pattern %.0f',i))
        for j = 1:nMovies
            [dataTracesFull, channelsFull] = NS_ReadPreprocessedDataNoArtifacts(DataPath,i,j);
            electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
            surroundingChannels = electrodeMap.getAdjacentsTo(i, 1)';
            dataTraces{i,j} = dataTracesFull(:,surroundingChannels,:);
            channels{i,j} = channelsFull(surroundingChannels);
        end
    end
end


%% clustering data

clusterIndex = cell(nPatterns, nMovies);
PCAScore = cell(nPatterns, nMovies);

for i = 1:nPatterns
    if i == 9 || i == 25 || i == 57
    else
        disp(sprintf('clustering pattern %.0f',i))
        for j = 1:nMovies
            nTraces = size(dataTraces{i,j}, 1);
            nElectrodes = size(dataTraces{i,j}, 2);
            nSamples = size(dataTraces{i,j}, 3);
            
            % concatenates traces on different electrodes for each pulse
            prinCompArray = zeros(nTraces, nElectrodes*nSamples);
            for k = 1:nElectrodes
                prinCompArray(:, (k-1)*nSamples + 1 : k*nSamples) = dataTraces{i,j}(:,k,:);
            end
            
            [PCACoef, PCAScore{i,j}] = princomp(prinCompArray,'econ');
            clusterIndex{i,j} = autoClusterImm(PCAScore{i,j}, verbose, iterations, DIMS, nBurnIn, a_0, b_0, mu_0, k_0, lambda_0, v_0);
        end
    end
end

fileName = sprintf('/Analysis.noindex/Lauren/2008-08-26-0/data006_clusters/clusterIndex_i%d_d%d_nB%d_a%.1f_b%.1f_k%.1f_v%.1f.mat ',...
    iterations, DIMS, nBurnIn, a_0, b_0, k_0, v_0);
save(fileName, 'clusterIndex')


%% plotting results

for i = 1:nPatterns
    if i == 9 || i == 25 || i == 57
    else
        for j = 1:nMovies
            figure
            hold on

            % plots traces as different colors representing the cluster they belong to.
            clusterColor = hsv(max(clusterIndex{i,j}));
            for k = 1:length(clusterIndex{i,j})
                plot(PCAScore{i,j}(k,1), PCAScore{i,j}(k,2), '.', 'MarkerFaceColor', clusterColor(clusterIndex{i,j}(k),:), 'MarkerEdgeColor', clusterColor(clusterIndex{i,j}(k),:), 'MarkerSize', 20)
            end

            hold off
        end
    end
end
figure
hold on
