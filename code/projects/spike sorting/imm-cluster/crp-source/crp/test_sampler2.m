% Copyright October, 2006, Brown University, Providence, RI. 
% All Rights Reserved 

% Permission to use, copy, modify, and distribute this software and its
% documentation for any purpose other than its incorporation into a commercial
% product is hereby granted without fee, provided that the above copyright
% notice appear in all copies and that both that copyright notice and this
% permission notice appear in supporting documentation, and that the name of
% Brown University not be used in advertising or publicity pertaining to
% distribution of the software without specific, written prior permission. 

% BROWN UNIVERSITY DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
% INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY
% PARTICULAR PURPOSE. IN NO EVENT SHALL BROWN UNIVERSITY BE LIABLE FOR ANY
% SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
% RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
% CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
% CONNECTION WITH THE USE.

% Author: Frank Wood fwood@cs.brown.edu.


addpath('/var/automount/netapp/snle/home/gauthier/Desktop/PCA testing/final version/crp/utilities')
addpath('/var/automount/netapp/snle/home/gauthier/Desktop/PCA testing/final version/crp/distributions')

% create Gaussian mixture model
g1 = gaussian([1 1]', .1*eye(2));
g2 = gaussian([-1 -1]', .1*eye(2));
g3 = gaussian([1 -1]', .1*eye(2));
g4 = gaussian([-1 1]', .1*eye(2));
gmm = gaussian_mixture_model(1/4,g1,1/4,g2,1/4,g3,1/4,g4);

% create training data by sampling from the Gaussian mixture model
%N = 1000;
% [test_data, test_data_class] = sample(gmm,N);
% training_data = test_data(:,1:5:1000);
% training_data_class = test_data_class(:,1:5:1000);

test_data = dataset{dd}.pca{1}.spike_projections(:,1:3)';
size(test_data)
training_data = test_data(:,1:1000);
%training_data_class = test_data_class(:,1:1000);


%class_label = unique(training_data_class);
size(training_data)
%training_data_class
%training_data
% show the training data
%plot_mixture(training_data,training_data_class)

points = 1:1000;
%points = 1000:2000;

spikeArray = training_data(:,points);
%spikeArray = spikes(points,1:5);


% run the CRP sampler to generate the posterior distribution over model 
% parameters
%[class_id, mean_record, covariance_record, K_record, lP_record,
%alpha_record, best] = sampler(spikeArray', 1000);
[class_id, mean_record, covariance_record, K_record, lP_record, alpha_record, best] = sampler(spikeArray, 200);
[class_id2, mean_record2, covariance_record2, K_record2, lP_record2, alpha_record2, best2] = sampler(spikeArray, 200);



%[class_id2, mean_record2, covariance_record2, K_record2, lP_record2, alpha_record2, best2] = sampler(spikeArray', 1000);
%[class_id3, mean_record3, covariance_record3, K_record3, lP_record3, alpha_record3, best3] = sampler(spikeArray', 1000);
figure(9); clf;  hold on; plot(K_record, 'r.');plot(K_record2, 'g.');%plot(K_record3, 'b.');
figure(10); clf; hold on; plot(lP_record,'r.');plot(lP_record2,'g.');%plot(lP_record3,'b.');

%plot map estimate ftw
%lP_best = [lP_record(best), lP_record2(best2)];%, lP_record3(best3)];

lP_best = [lP_record(best) lP_record2(best)];

%should always pick maximum, since -log likelihood value closest to zero is best
switch 3
    case 1;
        sorted = sort(lP_best);
        bestIteration = find(sorted(round(length(sorted/2))));
    case 2;
        bestIteration = find(min(lP_best));
    case 3;
        bestIteration = find(max(lP_best));
end;


figure(2); 
switch bestIteration
    case 1;
        bestIndex = find(lP_record == lP_best(bestIteration));
        %plot_mixture(spikes(points,1:2)', class_id(:,bestIndex));
        %plot_mixture(test_data(1:2,:), test_class(:,:));
    case 2; 
        bestIndex = find(lP_record2 == lP_best(bestIteration));
        %plot_mixture(spikes(points,1:2)', class_id(:,bestIndex2));
    case 3; 
        bestIndex = find(lP_record3 == lP_best(bestIteration));
        %plot_mixture(spikes(points,1:2)', class_id(:,bestIndex3));
end

figure(3);
switch bestIteration
    case 1;
        bestIndex = find(lP_record == lP_best(bestIteration));
        %plot_mixture(test_data(1:3,:), test_class(:,:));
        %plot_mixture(spikes(points,1:3)', class_id(:,bestIndex));
    case 2; 
        bestIndex = find(lP_record2 == max(lP_best(bestIteration)));
        %plot_mixture(spikes(points,1:3)', class_id(:,bestIndex2));
    case 3; 
        bestIndex = find(lP_record3 == max(lP_best(bestIteration)));
        %plot_mixture(spikes(points,1:3)', class_id(:,bestIndex3));
end

classes = zeros(size(training_data(:,1)),1);


mu_0 = zeros(size(training_data(:,1)));
k_0 = 1;
v_0 = size(training_data(:,1),1);
lambda_0 = eye(size(training_data(:,1),1))*.3;

lpa = zeros(K_record(best),1);

training_data = test_data;

nPoints = size(training_data,2);



%for each point
for point=1:nPoints
    %for each class
    for k=1:size(covariance_record{bestIndex},3)
        
        lp = lp_crp(k,alpha_record(bestIndex));
        
        sigma = covariance_record{bestIndex}(:,:,k);
        mu = mean_record{bestIndex}(:,k);

        lpa(k) = lp + fvnlp(training_data(:,point),mu, sigma);

    end

    %best cluster
    %lpa
    classes(point) = find(max(lpa) == lpa);
    lpa = zeros(K_record(best),1);
    lpa(lpa == 0) = NaN;
end

figure(555); plot_mixture(training_data(1:2,1:5:5000), classes(1:5:5000));  

%figure(666); hist(training_data_class);
%figure(999); hist(classes(1:nPoints));
