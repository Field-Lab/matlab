function positions = loadElecPositions512()
% Function that loads a matrix containing the electrode numbers for the 512-electrode MEA
currentFunctionPath = mfilename('fullpath');
temp = load(fullfile(currentFunctionPath, '../../resources/arrayPositions512.mat'));
positions = temp.positions;
end