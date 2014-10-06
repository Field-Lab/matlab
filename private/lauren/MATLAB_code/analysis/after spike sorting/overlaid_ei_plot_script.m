
pathToEi = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data001-lh/data001-lh.ei';
midgetIDs = [18 92 91 167 182 256 332 348 379 395 512 541 616 646 676 677 753 766 782 811 858 886 872 931];
nMidgets = length(midgetIDs);



[xCoords yCoords] = getElectrodeCoords61();


% get neuron's average spike waveforms from VISION analysis
eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);
eiAmps = cell(nMidgets, 1);
for i = 1:nMidgets
    ei = eiFile.getImage(midgetIDs(i)); %gets ei data for neuron,  as a 3D array: 2 x nElectrodes x nSamples
    eiAmps{i} = zeros(64, 1);
    for j = 1:64
        if ~(j==9||j==25||j==57)
            eiAmps{i}(j) = max(max(ei(1,j+1,:)));
        end
    end
    %eiAmpsNormalized = eiAmps/max(eiAmps);
end
clear eiFile

%removing unused electrodes

xCoords(57) = [];
xCoords(25) = [];
xCoords(9)  = [];
xCoords(4)  = [];
yCoords(57) = [];
yCoords(25) = [];
yCoords(9)  = [];
yCoords(4)  = [];
for i = 1:nMidgets
    eiAmps{i}(57) = [];
    eiAmps{i}(25) = [];
    eiAmps{i}(9) = [];
    eiAmps{i}(4) = [];  %electrode 4 is dead on first stim64 system
end

elecNos = [1:3 5:8 10:24 26:56 58:64];


%% contour preparation

eiAmpsGrid = zeros(9,9);
electrodeGrid = ...
    [0  0  40 42 45 48 51 0  0;...
     0  37 39 43 46 50 52 0  0;
     0  34 35 38 44 49 53 55 0;
     31 32 33 36 47 54 56 59 0;
     30 29 28 26 41 58 60 61 62;
     27 24 22 15 4  1  64 63 0;
     0  23 21 17 12 6  3  2  0;
     0  20 18 14 11 7  5  0  0;
     0  0  19 16 13 10 8  0  0];
 
 for i = 1:9
     for j = 1:9
         if electrodeOrder(i,j)
             eiAmpsGrid(i,j) = eiAmps(elecNos == electrodeGrid(i,j));
         end
     end
 end





figure
hold on
contour(xCoords, yCoords, eiAmps{1})


% figure
% hold on
% 
% for i = 1:64
%     if ~(i==9||i==25||i==57)
%         if eiAmpsNormalized(i)>.02
%             plot(xCoords(i), yCoords(i), 'ok','MarkerSize', ceil(eiAmpsNormalized(i)*20), 'MarkerFaceColor', 'k')
%         end
%     end
% end
% 
% plot(xCoords(pElec), yCoords(pElec), 'o', 'MarkerFaceColor', primColor, 'MarkerEdgeColor', primColor)
% plot(xCoords(sElec), yCoords(sElec), 'o', 'MarkerFaceColor', secColor, 'MarkerEdgeColor', secColor)
% %plot(xCoords(sElec)+0.4, yCoords(sElec), 'o', 'MarkerFaceColor', [1 0 1], 'MarkerEdgeColor', [1 0 1])


set(gca, 'XLim', [-10.4633 10.4633], 'YLim', [-8 8])
axis off
hold off

set(gca, 'XLim', [-10.4633 10.4633], 'YLim', [-8 8])