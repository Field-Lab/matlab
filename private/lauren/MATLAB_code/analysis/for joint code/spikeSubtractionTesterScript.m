clear all


addpath(genpath('/snle/home/lhruby/MATLAB_code/analysis/from Pawel/matlab_cp_2008_11_19/biblioteki/'))

details = getDatasetDetails(4);


nMovies = length(details.movieNumbers);
nPatterns = length(details.patternNumbers);


%% linkage-based artifact estimation

cd(details.savePath)
mkdir 'linkage_based'
cd 'linkage_based'


latencies = cell(length(details.patternNumbers), length(details.movieNumbers)); %pattern number x movie number
latenciesInitial = cell(length(details.patternNumbers), length(details.movieNumbers));

for i = 1:nMovies
    for j = 1:nPatterns
        
        [latencies{j,i} latenciesInitial{j,i}] = templateMatchClustering(details.dataPath, details.patternNumbers(j),...
            details.movieNumbers(i), details.pathToEi, details.neuronID, details.centerChannel, 'modelType', 'linkage');

    end
end

save latencies.mat latencies latenciesInitial

%% TTX arifact subtraction

if isfield(details, 'artifactMovieNumbers')

    cd(details.savePath)
    mkdir 'ttx_subtraction'
    cd 'ttx_subtraction'

    latencies = cell(length(details.patternNumbers), length(details.movieNumbers)); %pattern number x movie number
    latenciesInitial = cell(length(details.patternNumbers), length(details.movieNumbers));

    for i = 1:nMovies
        for j = 1:nPatterns
            
            ttx.path = details.artifactPath;
            ttx.patternNumber = details.artifactPatternNumbers(j);
            ttx.movieNumber = details.artifactMovieNumbers(i);
            
            [latencies{j,i} latenciesInitial{j,i}] = templateMatchClustering(details.dataPath, details.patternNumbers(j),...
                details.movieNumbers(i), details.pathToEi, details.neuronID, details.centerChannel, ttx, 'modelType', 'ttx');
            
            
%             [latencies{j,i} latenciesInitial{j,i}] = templateMatchClustering(details.dataPath, details.patternNumbers(j),...
%                 details.movieNumbers(i), details.pathToEi, details.neuronID, details.centerChannel, 0, details.artifactPath,...
%                 details.artifactPatternNumbers(j), details.artifactMovieNumbers(i));
        end
    end

    save latencies.mat latencies latenciesInitial

end

%% double exponential

cd(details.savePath)
mkdir 'sum_of_exponentials'
cd 'sum_of_exponentials'

latencies = cell(length(details.patternNumbers), length(details.movieNumbers)); %pattern number x movie number
latenciesInitial = cell(length(details.patternNumbers), length(details.movieNumbers));


for i = 1:nMovies
    for j = 1:nPatterns
        
        [latencies{j,i} latenciesInitial{j,i}] = templateMatchClustering(details.dataPath, details.patternNumbers(j),...
            details.movieNumbers(i), details.pathToEi, details.neuronID, details.centerChannel, 'modelType', 'sumExp');
        
        
%         [latencies{j,i} latenciesInitial{j,i}] = templateMatchClustering(details.dataPath, details.patternNumbers(j),...
%             details.movieNumbers(i), details.pathToEi, details.neuronID, details.centerChannel);
    end
end

save latencies.mat latencies latenciesInitial

%% previous artifact

cd(details.savePath)
mkdir 'previous_artifact'
cd 'previous_artifact'

latencies = cell(length(details.patternNumbers), length(details.movieNumbers)); %pattern number x movie number
latenciesInitial = cell(length(details.patternNumbers), length(details.movieNumbers));


prevArtifactVector = [];
for i = 1:nMovies
    for j = 1:nPatterns
        [latencies{j,i} latenciesInitial{j,i} prevArtifactVector] = templateMatchClustering(details.dataPath, details.patternNumbers(j),...
            details.movieNumbers(i), details.pathToEi, details.neuronID, details.centerChannel, 'modelType', 'prevArtifact', 'prevArtifact', prevArtifactVector);
        
%         [latencies{j,i} latenciesInitial{j,i} prevArtifact] = templateMatchClustering(details.dataPath, details.patternNumbers(j),...
%             details.movieNumbers(i), details.pathToEi, details.neuronID, details.centerChannel, prevArtifact);
    end
end

save latencies.mat latencies latenciesInitial

