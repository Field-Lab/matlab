% For each neuron, plot the sensitivity map. 
eiFilePath = '/Volumes/Analysis/2015-11-09-3/data000/data000.ei'; 
fPath = '/Volumes/Analysis/2015-11-09-3/data001-autosort/';
[eiM,neuronIdList] = convertEiFileToMatrix(eiFilePath);
figure(101); 
for n = 1:length(neuronIdList)
    cla;
    neuronId = neuronIdList(n);
    eiAmps = max(eiM(:,:,n),[],2) - min(eiM(:,:,n),[],2); 
    eiContour_wLinFit(eiAmps,'figureNum',101); 
    title(['Neuron ' num2str(neuronId)]); 
    [allVals, actProb, actAmp] = genActThreshSpatialMaps_auto(fPath,neuronId);
    idx = find(allVals); 
    scatter(xc(idx),yc(idx),100,allVals(idx),'filled'); axis image; colorbar;
    pause;
end