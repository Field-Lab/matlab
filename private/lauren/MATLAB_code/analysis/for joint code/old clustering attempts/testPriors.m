addpath('/snle/lab/Development/RRS/imm-cluster/crp-source/utilities');
addpath('/snle/lab/Development/RRS/imm-cluster/crp-source/distributions');
d = pwd;

% Change directory containing IMM/CRP code
cd('/snle/lab/Development/RRS/imm-cluster/crp-source/crp/')

% Parameters for the gamma function (sampled from when computing 
% alpha in the chinese restaurant process sampler)
a_0      = 1;
b_0      = 1;

% Priors on cluster means
mu_0     = zeros(nDimensions,1);
k_0      = .1;

% Priors on cluster covariances
lambda_0 = eye(nDimensions)* 1000;
v_0      = 10;

showPlots = true;   %whether or not to plot generated data at each iteration of loop


%% find mean and variance in number of clusters for different values of chosen parameter

meanClusters = zeros(1,100);
varianceClusters = zeros(1,100);

for k = 1:100
    k
    iteratedClusters = zeros(1,50);
    close all
    b_0 = 0.01*k;
    
    for i = 1:50

        n = 100; %number of data points to be generated
        [datapoints labels] = sample_igmm_prior(n,a_0,b_0,mu_0,lambda_0,k_0,v_0);

        clusters = length(unique(labels));
        sampleColors = hsv(clusters);
        if showPlots
            figure(1)
            plot(datapoints(:,1), datapoints(:,2), '.')
            title(sprintf('sampled from priors; number of clusters = %d', clusters))
            iteratedClusters(i) = clusters;
            pause(0.1)
        end
    end
    meanClusters(k) = mean(iteratedClusters);
    varianceClusters(k) = var(iteratedClusters);
end

figure
hold on
subplot(2,1,1)
plot(0.01:0.01:1, meanClusters)
subplot(2,1,2)
title('mean number of clusters')
plot(0.01:0.01:1, varianceClusters)
title('variance in number of clusters')
xlabel('value of a_0')
hold off



% change back to original directory
cd (d);