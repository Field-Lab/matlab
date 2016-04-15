%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stuff to change for each piece

piece = '2016-02-17-6';
Analysis_Path = ['/Volumes/Analysis/' piece '/mVision/'];
fig_save = ['/Users/Nora/Desktop/Fig_Output/' piece '/mVision/'];

% dataruns
class = 'data000';
% split = 'data005';
reg{1} = 'data003';

% Maskin
%%{
fig_save = [fig_save 'A'];
masking = 'data002';
cells_orig = [6023, 1052, 3421, 7522, 4653, 2221];
%cells_reg{1} = cells_orig;
cells_masking = cells_orig;
cells_reg{1} = [6017;1052;3423;7518;4652;2222];
cells_masking = [6139; 1053; 3423; 7521; 4653; 2222];
n_masks = 7;
sigmas = [2 2 4 4 5 5 6 6 4 4 5 5 6 6];
% subgroup 1
cell_idx{1} = 1:3;
mask_conditions{1} = [1 3 5 7];
comp_conditions{1} = [2 4 6 8];
% subgroup 1
cell_idx{2} = 4:6;
mask_conditions{2} = [1 9 11 13];
comp_conditions{2} = [2 10 12 14];
%}
%{
masking = 'data004';
fig_save = [fig_save 'B'];
cells_orig = [5929, 1113, 3668 7414, 4592, 2342];
cells_reg{1} = [5929;1112;3671;7276;4699;2342];
cells_masking = [6034;1112;3666;7280;4715;2361];;
n_masks = 7;
sigmas = [2 2 4 4 5 5 6 6 4 4 5 5 6 6];
% subgroup 1
cell_idx{1} = 1:3;
mask_conditions{1} = [1 3 5 7];
comp_conditions{1} = [2 4 6 8];
% subgroup 1
cell_idx{2} = 4:6;
mask_conditions{2} = [1 9 11 13];
comp_conditions{2} = [2 10 12 14];
%}
