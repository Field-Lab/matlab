patternsNeu1699=[477:515 556:561];

for p=1:length(patternsNeu1699)
    
[thresholdHum thresholdAlg] = fitToErfOutputAndHuman(Outputs10Jun(patternsNeu1699(p)));
patternNos(p) = Outputs10Jun(patternsNeu1699(p)).stimInfo.patternNo;
thresHumNeu1699(p) = thresholdHum;
thesholdAlgNeu1699(p) = thresholdAlg;
end

allValsH = zeros(512,1); 
allValsH(patternNos) = thresHumNeu1699; 
allValsA = zeros(512,1); 
allValsA(patternNos) = thesholdAlgNeu1699; 
humanVals = thresHumNeu1699;
algorithmVals = thesholdAlgNeu1699;
algorithmVals(algorithmVals>4) = 0; 
algorithmVals(algorithmVals<0) = 0; 
humanVals(humanVals>4) = 0; 
humanVals(humanVals<0) = 0; 
figure; line(0:4,0:4,'Color','k'); hold on; 
scatter(humanVals, algorithmVals,100,'filled');
xlabel('thresholds (\muA) using human analysis'); 
ylabel('thresholds (\muA) using GM algorithm'); 
title('50% activation thresholds for n1699'); 
%% Plot the EI contour with a linear fit + electrodes
[xc,yc] = getElectrodeCoords512(); 

temp = load('/Volumes/Analysis/2015-04-14-0/data001/elecResp_n1699_p410.mat');
ei = temp.elecResp.cells.mainEI;
neuronId = temp.elecResp.cells.main;
clear temp; 
eiAmps = max(ei,[],2) - min(ei,[],2);
eiContour_wLinFit(eiAmps,'linFitThresh',5)
title(sprintf('neuron %d',neuronId)); 

% Show stimulating electrodes.
trycmap = ([linspace(70,156,64); linspace(130,49,64); linspace(180,8,64)]/255)'; 
cmap = flipud(hot); 
cmap = load('/Users/grosberg/matlab/dataset_specific/redtealcmap.mat'); 
allValsA(allValsA>4) = 0; 
allValsA(allValsA<0) = 0; 
idx = find(allValsA); 
hold on; scatter(xc(idx),yc(idx),150, allValsA(idx),'filled'); 
colormap(cmap.m); caxis([0.5 3.5]); colorbar; 

eiContour_wLinFit(eiAmps,'linFitThresh',5)
title(sprintf('neuron %d',neuronId)); 

% Show stimulating electrodes.
 
cmap = load('/Users/grosberg/matlab/dataset_specific/redtealcmap.mat'); 
allValsH(allValsH>4) = 0; 
allValsH(allValsH<0) = 0; 
idx = find(allValsH); 
hold on; scatter(xc(idx),yc(idx),150, allValsH(idx),'filled'); 
colormap(cmap.m); caxis([0.5 3.5]); colorbar; 