elecStimPath = '/Volumes/Analysis/2014-08-20-1/data003/'; %'/Volumes/Analysis/2012-09-24-3/data008/';

patternNos = [65 69 73 77 81 89 259 263 264 267 271 272 275 279 ...
    280 283 284 288 291 292 296 300 304 308 316];
display = false;
patternNos = [65:128 257:320];
bundleMaxes = zeros(368, 1);
for patternNo = patternNos
    bundleMean = getBundleVoltagesAStar(patternNo, display, 0);
    bundleMaxes(patternNo) = max(abs(bundleMean(:,1)));
end
disp('finished analyzing patterns');
electrodesByPattern = getElecsFromPatternFiles([elecStimPath 'pattern_files/']);

electrodeNumberingScheme = zeros(size(electrodesByPattern, 1), 1);
electrodeNumberingScheme(patternNos) = electrodesByPattern(patternNos);
axonVoltageByElectrode = zeros(1, 512);

%axonVoltageByElectrode(electrodeNumberingScheme(patternNos)) = bundleMaxes(patternNos);

axonVoltageByElectrode(electrodeNumberingScheme(patternNos)) = bundleMaxes(patternNos);
%axonVoltageByElectrode(electrodeNumberingScheme(biElecPatterns1)) = biElec1Maxes(biElecPatterns1);
%axonVoltageByElectrode(electrodeNumberingScheme(biElecPatterns2)) = biElec2Maxes(biElecPatterns2);


[~, electrodeArray] = ei2matrix(axonVoltageByElectrode);
cmap = colormap(flipud(gray));
%cmap(1,:) = [0,0,0]

electrodeArray = electrodeArray / 16.05148356;

figure; imagesc(electrodeArray); colorbar; colormap(cmap); title('Approx. # of axons - bi-electrode 2'); axis image; set(gcf, 'InvertHardCopy', 'Off');


%%

for patternNo = patternNos
    plotBundleVoltage(patternNo, false);
    saveas(gcf, ['p' num2str(patternNo)], 'epsc');
end