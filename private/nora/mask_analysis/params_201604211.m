%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stuff to change for each piece
clear reg
piece = '2016-04-21-1';
Analysis_Path = ['/Volumes/Analysis/' piece '/mVision/'];
fig_save = ['/Volumes/Lab/Users/Nora/Fig_Output/' piece '/'];

% dataruns
class = ['/Volumes/Analysis/' piece '/streamed/data015/data015'];
reg{1} = 'data018';

% Maskin
%%{
masking = 'data017';
cells_orig = [811 3169 6047 2281 2341 5401];
cells_reg{1} = [812;3170;6054;2286;2347;4624];
cells_masking = [813;3287;6048;NaN;2341;4501];
n_masks = 7;
sigmas = [2 2 4 4 6 6 8 8 10 10 8 8 10 10];
% subgroup 1
cell_idx{1} = 1:3;
cell_type{1} = 'on';
mask_conditions{1} = [1 3 5 7 9];
comp_conditions{1} = [2 4 6 8 10];
% subgroup 1
cell_idx{2} = 4:6;
cell_type{2} = 'off';
mask_conditions{2} = [1 3 5 11 13];
comp_conditions{2} = [2 4 6 12 14];
%}

%{
% Maskin2
% WARNING THE TIMING IS MESSED UP ON THIS ONE!!!
fig_save = [fig_save 'B'];
masking = 'data023';
%cells_orig = [1369 2299 3394 4351 6046 7367];
cells_orig = [1369;2297;3393;4353;6048;7366];
cells_reg{1} = [1368;2296;3512;4352;6047;7366];
cells_masking = cells_orig;
n_masks = 7;
sigmas = [2 2 4 4 5 5 6 6 4 4 5 5 6 6];
% subgroup 1
cell_idx{1} = 2:2:6;
mask_conditions{1} = [1 3 5 7];
comp_conditions{1} = [2 4 6 8];
% subgroup 1
cell_idx{2} = 1:2:6;
mask_conditions{2} = [1 9 11 13];
comp_conditions{2} = [2 10 12 14];
%}

% for i = 1:n_groups
%     mask_order{i} = ['/Volumes/Data/' piece '/Visual/masks/' masking{i} '_mask_files.mat'];
% end