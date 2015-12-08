filePath = '/Volumes/Analysis/2015-10-06-3/data004/';
elecsUsed = getElecsFromPatternFiles([filePath 'pattern_files/']);
figure; imagesc(elecsUsed); 
% elecsOfInterest
% 199/200: p 13/14
% 207/208: p 29/30 p29 appears to be activating a single axon at movieIndex 31
% 214/215: p 43/44
% 222/223: p 59/60
% 231/232: p 77/78
% 239/240: p 93/94
% 307/315: p 229/230

patternNos2e = [13 14 29 30 43 44 59 60 77 78 93 94 229 230]; 
patternNos1e=[247 248 255 256 262 263 270 271 279 280 287 288 355 363]; 
patternNo = patternNos2e(14); 
playMovie512arrayAfterStimPattern_dots(filePath,patternNo,'movieIndex',10); 

thresholds2e = [1.1 1.41 2.01 1.81 1.04 1.2 1.1 1.2 1.1 1.41 1.71 1.41 1.51 1.71];
patternNo = patternNos1e(14); 
thresholds1e = [1.1 1.71 2.21 1.81 1.1 1.1 1.04 1.2 1.04 1.41 1.81 1.41 1.41 1.71];

figure; plot(thresholds1e,'-xb'); 
hold on; plot(thresholds2e,'-om'); legend('1 elec stim','2 elec stim');
title('Comparison of axon bundle activation thresholds for 1- and 2-electrode stimulation'); 
for p = 1:length(patternNos1e)
    patternNo = patternNos1e(p);
    ms = findMovieNos(filePath,patternNo);    
    [~, elecs1e(p), ~, ~, ~, ~] = getStimAmps(filePath, patternNo, ms(1));
end

for p = 1:length(patternNos2e)
    patternNo = patternNos2e(p);
    ms = findMovieNos(filePath,patternNo);    
    [amps, elecs, ~, ~, ~, ~] = getStimAmps(filePath, patternNo, ms(1));
    elecs2e(p)= elecs(find(amps<0));
end
elecs1e == elecs2e
figure; plot(elecs1e,thresholds1e,'-xb'); 
hold on; plot(elecs2e,thresholds2e,'-om'); legend('1 elec stim','2 elec stim');
title('Comparison of axon bundle activation thresholds for 1- and 2-electrode stimulation'); 
xlabel('stimulating electrode'); 