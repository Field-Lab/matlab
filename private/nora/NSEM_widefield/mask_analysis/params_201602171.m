%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stuff to change for each piece

piece = '2016-02-17-1';
Analysis_Path = ['/Volumes/Analysis/' piece '/mVision/'];
fig_save = ['/Volumes/Lab/Users/Nora/Fig_Output/' piece '/mVision/'];

% dataruns
class = ['/Volumes/Analysis/' piece '/streamed/data019/data019'];
split = 'data020';
reg{1} = 'data022';
long = 'data025';

% Maskin
%%{
masking = 'data021';
fig_save = [fig_save 'A'];
%cells_orig = [7172, 4805, 2419, 5972, 1051, 3316];
cells_orig = [7172, 4804, 2401, 5971, 1051, 3316];
cells_reg{1} = [7172;4805;2416;5972;1053;3437];
cells_masking = [7171;4863;2416;5957;1052;3318];
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