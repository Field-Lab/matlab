%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stuff to change for each piece

piece = '2016-02-17-6';
Analysis_Path = ['/Volumes/Analysis/' piece '/mVision/'];
fig_save = ['/Volumes/Lab/Users/Nora/Fig_Output/' piece '/mVision/'];

% dataruns
class =['/Volumes/Analysis/' piece '/streamed/data000/data000'];
% split = 'data005';
reg{1} = 'data003';

% Maskin
%%{
fig_save = [fig_save 'A'];
masking = 'data002';
cells_orig = [6023, 1052, 3421, 7522, 4653, 2221];

cells_reg{1} = [6017;1052;3423;7518;4652;2222];
cells_masking = [6139; 1053; 3423; 7521; 4653; 2222];
n_masks = 7;
sigmas = [2 2 4 4 5 5 6 6 4 4 5 5 6 6];
% subgroup 1
%off
cell_type{1} = 'off';
cell_idx{1} = 1:3;
mask_conditions{1} = [1 3 5 7];
comp_conditions{1} = [2 4 6 8];
% subgroup 1
%on
cell_idx{2} = 4:6;
cell_type{2} = 'on';
mask_conditions{2} = [1 9 11 13];
comp_conditions{2} = [2 10 12 14];
%}

