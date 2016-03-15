% Histogram of the bundle thresholds. 
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/hand_sorted_bundle_thresholds/axonBundleThresholds_byPattern_2012_09_24_3_data001.mat')
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/hand_sorted_bundle_thresholds/axonBundleThresholds_2015_10_06_3.mat')
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/hand_sorted_bundle_thresholds/axonBundleThresholds_byPattern_2015_04_09_2_data003.mat')
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/hand_sorted_bundle_thresholds/axonBundleThresholds_2015_10_06_6.mat')
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/hand_sorted_bundle_thresholds/axonBundleThresholds_2015_09_23_2.mat')


allthresholds = cat(1,axonBundleThresholds',axonBundleThresholds_2015_09_23_2,axonBundleThresholds_2015_10_6_3,bundleThresholds_2015_10_06_6,bundle_thresholds_2012_09_24); 
allthresholds(allthresholds==0) = []; 
allthresholds(isnan(allthresholds)) = []; 
figure; histogram(allthresholds/1.1,'binmethod','sturges'); 
xlabel('stimulation current amplitude (uA)'); 
ylabel('number of measurements');
title(sprintf('Largest current applied without activating in axon bundle in\n2012-09-24-3, 2015-04-09-2, 2015-09-23-2, 2015-10-06-3, 2015-10-06-6')); 
figure; 
h1 = histfit(allthresholds,13,'kernel');