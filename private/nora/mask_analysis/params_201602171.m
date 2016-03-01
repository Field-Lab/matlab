%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stuff to change for each piece

piece = '2016-02-17-1';

% dataruns
class = 'data019';
masking{1} = 'data021';
masking{2} = 'data023';
split = 'data020';
reg{1} = 'data022';

% Maskin
cells{1} = [7172, 4804, 2401, 5971, 1051, 3316];
n_masks{1} = 7;
sigmas{1} = [2 2 4 4 5 5 6 6 4 4 5 5 6 6];
% subgroup 1
cell_idx{1}{1} = 1:3;
mask_conditions{1}{1} = [1 3 5 7];
comp_conditions{1}{1} = [2 4 6 8];
% subgroup 1
cell_idx{1}{2} = 4:6;
mask_conditions{1}{2} = [1 9 11 13];
comp_conditions{1}{2} = [2 10 12 14];

% Maskin2
cells{2} = [7366, 4351, 2296, 6046, 1366, 3393];
n_masks{2} = 7;
sigmas{2} = [2 2 4 4 5 5 6 6 4 4 5 5 6 6];
% subgroup 1
cell_idx{2}{1} = 1:3;
mask_conditions{2}{1} = [1 3 5 7];
comp_conditions{2}{1} = [2 4 6 8];
% subgroup 1
cell_idx{2}{2} = 4:6;
mask_conditions{2}{2} = [1 9 11 13];
comp_conditions{2}{2} = [2 10 12 14];

% for i = 1:n_groups
%     mask_order{i} = ['/Volumes/Data/' piece '/Visual/masks/' masking{i} '_mask_files.mat'];
% end