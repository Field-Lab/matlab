%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stuff to change for each piece

piece = '2016-02-17-6';
Analysis_Path = ['/Volumes/Analysis/' piece '/mVision/'];
fig_save = ['/Volumes/Lab/Users/Nora/Fig_Output/' piece '/mVision/'];

% dataruns
class =['/Volumes/Analysis/' piece '/streamed/data000/data000'];
% split = 'data005';
reg{1} = 'data003';

%%{
masking = 'data004';
fig_save = [fig_save 'B'];
cells_orig = [5929, 1113, 3668 7414, 4592, 2342];
cells_reg{1} = [5929;1112;3671;7276;4699;2342];
cells_masking = [6034;1112;3666;7280;4715;2361];
n_masks = 7;
sigmas = [2 2 4 4 5 5 6 6 4 4 5 5 6 6];
% subgroup 1
cell_idx{1} = 1:3;
cell_type{1} = 'on';
mask_conditions{1} = [1 3 5 7];
comp_conditions{1} = [2 4 6 8];
% subgroup 1
cell_idx{2} = 4:6;
cell_type{2} = 'off';
mask_conditions{2} = [1 9 11 13];
comp_conditions{2} = [2 10 12 14];
%}