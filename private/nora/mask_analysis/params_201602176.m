%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stuff to change for each piece

piece = '2016-02-17-6';

% dataruns
class = 'data000';
masking{1} = 'data002';
masking{2} = 'data004';
split = 'data005';
reg{1} = 'data003';

% Maskin
cells{1} = [6023, 1052, 3421 7522, 4653, 2221];
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
cells{2} = [5929, 1113, 3668 7414, 4592, 2342];
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