function plotElecsDots(stimelec,elecs)
%%Given list of electrodes, plots em on array

positions = loadElecPositions512();
f = figure; set(f,'Position',[100 465 845 445]);
set(f,'Color','white');
filld = zeros(1,512); filld(elecs) = 1; filld(stimelec) = 2;
scatter(positions(:,1), positions(:,2), 350, filld, 'filled')
