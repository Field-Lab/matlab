i = 1;
estimates = [];
patterns = [1:218];
for i = patterns(1:end)
    
    bundleMeans = abs(getBundleVoltagesAStar('/Volumes/Analysis/2015-05-27-0/data002/', i, false));
    avg = mean(bundleMeans(1:5, 1));
    index = find(bundleMeans(:,1) > 20, 1, 'first');
    %index = find(bundleMeans(:,1) > 5 * avg, 1, 'first');
    disp(num2str(bundleMeans(index, 3)));
    bMean = bundleMeans(index, 3);
    if isempty(bMean)
        bMean = 0;
    end
    estimates = [estimates; bMean];
    
end
