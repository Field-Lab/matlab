% Histogram of the bundle thresholds. 
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/hand_sorted_bundle_thresholds/axonBundleThresholds_byPattern_2012_09_24_3_data001.mat')
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/hand_sorted_bundle_thresholds/axonBundleThresholds_2015_10_06_3.mat')
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/hand_sorted_bundle_thresholds/axonBundleThresholds_byPattern_2015_04_09_2_data003.mat')
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/hand_sorted_bundle_thresholds/axonBundleThresholds_2015_10_06_6.mat')
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/hand_sorted_bundle_thresholds/axonBundleThresholds_2015_09_23_2.mat')

axonBundleThresholds(axonBundleThresholds==0) = []; 
axonBundleThresholds(isnan(axonBundleThresholds)) = []; 
axonBundleThresholds_2015_09_23_2(axonBundleThresholds_2015_09_23_2==0) = []; 
axonBundleThresholds_2015_09_23_2(isnan(axonBundleThresholds_2015_09_23_2)) = []; 
axonBundleThresholds_2015_10_6_3(axonBundleThresholds_2015_10_6_3==0) = []; 
axonBundleThresholds_2015_10_6_3(isnan(axonBundleThresholds_2015_10_6_3)) = []; 
bundleThresholds_2015_10_06_6(bundleThresholds_2015_10_06_6==0) = []; 
bundleThresholds_2015_10_06_6(isnan(bundleThresholds_2015_10_06_6)) = []; 
bundle_thresholds_2012_09_24(bundle_thresholds_2012_09_24==0) = []; 
bundle_thresholds_2012_09_24(isnan(bundle_thresholds_2012_09_24)) = []; 

allthresholds = cat(1,axonBundleThresholds',axonBundleThresholds_2015_09_23_2,axonBundleThresholds_2015_10_6_3,bundleThresholds_2015_10_06_6,bundle_thresholds_2012_09_24); 
allthresholds(allthresholds==0) = []; 
allthresholds(isnan(allthresholds)) = []; 
figure; h=histogram(allthresholds/1.1,'binmethod','sturges'); 
xlabel('stimulation current amplitude (uA)'); 
ylabel('number of measurements');
title(sprintf('Largest current applied without activating in axon bundle in\n2012-09-24-3, 2015-04-09-2, 2015-09-23-2, 2015-10-06-3, 2015-10-06-6')); 
figure; 
h1 = histfit(allthresholds,13,'kernel');

names = [repmat('2015_04_09_2',length(axonBundleThresholds'),1);
    repmat('2015_09_23_2',length(axonBundleThresholds_2015_09_23_2),1);
    repmat('2015_10_06_3',length(axonBundleThresholds_2015_10_6_3),1);
    repmat('2015_10_06_6',length(bundleThresholds_2015_10_06_6),1);
    repmat('2012_09_24_3',length(bundle_thresholds_2012_09_24),1) ];

figure; boxplot(allthresholds/1.1,names); 
ylabel('stimulation current amplitude (\muA)'); 
title('safe zone'); 

% Try stacked histogram. 
[n,edges,bin] = histcounts(allthresholds/1.1,'binmethod','sturges');
% edges = h.BinEdges; 
n1 = histcounts(axonBundleThresholds'/1.1,edges);
n2 = histcounts(axonBundleThresholds_2015_09_23_2/1.1,edges);
n3 = histcounts(axonBundleThresholds_2015_10_6_3/1.1,edges);
n4 = histcounts(bundleThresholds_2015_10_06_6/1.1,edges);
n5 = histcounts(bundle_thresholds_2012_09_24/1.1,edges);

figure; bar(edges(1:end-1)+diff(edges)/2,[n1; n2; n3; n4; n5]',1,'stacked'); 
xlabel('stimulation current amplitude (uA)'); 
ylabel('number of measurements');legend('retina 1','retina 2','retina 3','retina 4', 'retina 5'); 


%% Calculate the mean, median, and std deviattion of the results

figure; 
datasets={'2015-04-09-2','2015-09-23-2','2015-10-06-3','2015-10-06-6','2012-09-24-3'}; 
for d = 1:length(datasets);
    dataset = datasets{d};
    switch dataset
        case '2015-04-09-2'
            test_data = axonBundleThresholds';
            remaining = cat(1,axonBundleThresholds_2015_09_23_2,axonBundleThresholds_2015_10_6_3,...
                bundleThresholds_2015_10_06_6,bundle_thresholds_2012_09_24);
        case '2015-09-23-2'
            test_data = axonBundleThresholds_2015_09_23_2;
            remaining = cat(1,axonBundleThresholds',axonBundleThresholds_2015_10_6_3,...
                bundleThresholds_2015_10_06_6,bundle_thresholds_2012_09_24);     
        case '2012-09-24-3'
            test_data = bundle_thresholds_2012_09_24;
            remaining = cat(1,axonBundleThresholds',axonBundleThresholds_2015_09_23_2,...
                axonBundleThresholds_2015_10_6_3,bundleThresholds_2015_10_06_6);
        case '2015-10-06-6'
            test_data = bundleThresholds_2015_10_06_6;
            remaining = cat(1,axonBundleThresholds',axonBundleThresholds_2015_09_23_2,...
                axonBundleThresholds_2015_10_6_3,bundle_thresholds_2012_09_24);
        case '2015-10-06-3'
            test_data =axonBundleThresholds_2015_10_6_3;
            remaining =cat(1,axonBundleThresholds',axonBundleThresholds_2015_09_23_2,...
                bundleThresholds_2015_10_06_6,bundle_thresholds_2012_09_24);
    end
   subplot(2,3,d); boxplot([test_data; remaining],[repmat(dataset,length(test_data),1);
        repmat('allOtherData',length(remaining),1) ]);
    [h,p,ci,stats] = ttest2(test_data,remaining);
    if h == 0
        title('means are the same')
    elseif h==1
        title('means not the same')
    end
end
%% Do resampling analysis. 
allT = allthresholds; 
allT(allT>3.7) = []; 
y = datasample(allthresholds,size(allthresholds,1),'Replace',true); 
[bootstat,bootsam] = bootstrp(1000,@mean,allthresholds);
[fi,xi] = ksdensity(bootstat);
se = std(bootstat); 
figure; subplot(1,2,1); 
histogram(bootstat);
xlabel('mean of all thresholds \muA'); ylabel('number of bootstrap samples'); 
subplot(1,2,2); plot(xi,fi); hold on; plot(xi(find(fi==max(fi))),max(fi),'ro'); 
xlabel('mean of all thresholds \muA'); ylabel('kernel density of bootstrapped samples'); 
fprintf('the max kernel density of 1000 bootstrapped means is at %0.2f uA with %0.3f standard error\n',...
    xi(find(fi==max(fi))),se);

[bootstat,bootsam] = bootstrp(1000,@median,allthresholds);
[fi,xi] = ksdensity(bootstat);
se = std(bootstat); 
figure; subplot(1,2,1); 
histogram(bootstat);
xlabel('median of all thresholds \muA'); ylabel('number of bootstrap samples'); 
subplot(1,2,2); plot(xi,fi);
xlabel('median of all thresholds \muA'); ylabel('kernel density of bootstrapped samples'); 

fprintf('the max kernel density of 1000 bootstrapped medians is at %0.2f uA with %0.3f standard error\n',...
    xi(find(fi==max(fi))),se);
%% Test the difference in means. 
nReps = 10000;
aa = .05;        %alpha value
datasets={'2015-04-09-2','2015-09-23-2','2015-10-06-3','2015-10-06-6','2012-09-24-3'}; 
for d = 1:length(datasets);
    dataset = datasets{d};
    switch dataset
        case '2015-04-09-2'
            x1 = axonBundleThresholds';
        case '2015-09-23-2'
            x1 = axonBundleThresholds_2015_09_23_2;
        case '2012-09-24-3'
            x1 = bundle_thresholds_2012_09_24;
        case '2015-10-06-6'
            x1 = bundleThresholds_2015_10_06_6;
        case '2015-10-06-3'
            x1 =axonBundleThresholds_2015_10_6_3;
    end
    
    x2 = allthresholds;
    n1 = size(x1,1);            %sample size 1
    n2 = size(x2,1);            %sample size 2
    
    % define the test statistic as the difference between means
    myStatistic = @(x1,x2) mean(x1)-mean(x2);
    
    sampStat = myStatistic(x1,x2);
    bootstrapStat = zeros(nReps,1);
    for i=1:nReps
        sampX1 = x1(ceil(rand(n1,1)*n1));
        sampX2 = x2(ceil(rand(n2,1)*n2));
        bootstrapStat(i) = myStatistic(sampX1,sampX2);
    end
    % Calculate the confidence interval (I could make a function out of this...)
    CI = prctile(bootstrapStat,[100*aa/2,100*(1-aa/2)]);
    
    %Hypothesis test: Does the confidence interval cover zero?
    H = CI(1)>0 | CI(2)<0;
    figure;
    xx = min(bootstrapStat):.01:max(bootstrapStat);
    hist(bootstrapStat,xx);
    hold on
    ylim = get(gca,'YLim');
    h1=plot(sampStat*[1,1],ylim,'y-','LineWidth',2);
    h2=plot(CI(1)*[1,1],ylim,'r-','LineWidth',2);
    plot(CI(2)*[1,1],ylim,'r-','LineWidth',2);
    h3=plot([0,0],ylim,'b-','LineWidth',2);
    xlabel('Difference between means');
    
    decision = {'Fail to reject H0','Reject H0'};
    title([dataset ' ' decision(H+1)]);
    legend([h1,h2,h3],{'Sample mean',sprintf('%2.0f%% CI',100*aa),'H0 mean'},'Location','NorthWest');
end