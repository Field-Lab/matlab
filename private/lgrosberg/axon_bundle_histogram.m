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