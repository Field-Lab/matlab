temp = load('~/git_code/matlab/private/lauren/MATLAB_code/analysis/dataset-specific/axonBundleThresholds_byPattern_2012_09_24_3_data008.mat');
axonBundleThresholds = temp.axonBundleThresholds_byPattern_2012_09_24_3_data008;

i = 1;
estimates = [];
for i = var(18:end)
    
    bundleMeans = abs(getBundleVoltagesAStar('/Volumes/Analysis/2015-04-09-2/data003/', i, false));
    avg = mean(bundleMeans(1:5, 1));
    index = find(bundleMeans(:,1) > 20, 1, 'first');
    %index = find(bundleMeans(:,1) > 5 * avg, 1, 'first');
    disp(num2str(bundleMeans(index, 3)));
    
    estimates = [estimates; bundleMeans(index, 3)];
    
end