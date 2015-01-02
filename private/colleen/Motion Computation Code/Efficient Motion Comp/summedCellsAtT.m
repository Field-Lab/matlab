%This function computes the equation N = sum_t(sum_i(r_i(t-\Delta t))^2) -
%sum_t(sum_i(r_i(t+\Delta t))^2) and returns each term in that equation
%separately.
% --------
%Inputs:
% time:
% flt_rsp2_shiftedRight
% flt_rsp2_shiftedLeft
% indices1
% spks_2

%Outputs:
% summedCellsR : motion signal from rightward shift
% summedCellsL : motion signal from leftward shift


function [summedCellsR, summedCellsL] = summedCellsAtT(time, flt_rsp2_shiftedRight, flt_rsp2_shiftedLeft, indices1, spks_2)
% make time compatible with cell functions
time = num2cell((repmat(time, size(indices1,2),1))) ;
time_cell = cell(size(indices1'));
time_cell = cellfun(@(x,t) t, time_cell, time, 'UniformOutput', false);
% Find cells that didn't spike
emptyCells = cellfun(@isempty,spks_2);

% Gaussian filtering of every cell 
everyCellAtTRight = cell2mat(cellfun(@(x,t) x(t), flt_rsp2_shiftedRight,time_cell, 'UniformOutput', false));
everyCellAtTLeft = cell2mat(cellfun(@(x,t) x(t), flt_rsp2_shiftedLeft,time_cell, 'UniformOutput', false));

% Adjust for cells that never fired
everyCellAtTLeft(emptyCells) = 0;
everyCellAtTRight(emptyCells) = 0;

% Sum to compute output
summedCellsR = sum(everyCellAtTRight);
summedCellsL = sum(everyCellAtTLeft);
