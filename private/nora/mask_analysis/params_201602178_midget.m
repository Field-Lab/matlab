%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stuff to change for each piece

piece = '2016-02-17-8';

class = 'data006';
masking{1} = 'data008';
split = 'data001';
reg{1} = 'data010';

% Midgets
cells{1} = [1157, 1561, 2388, 3031, 3977, 4867, 5451, 5956, 6458, 7265];
n_masks{1} = 2;
sigmas{1} = [2 2 4 4];
% subgroup 1
cell_idx{1}{1} = 1:10;
mask_conditions{1}{1} = [1 3];
comp_conditions{1}{1} = [2 4];

% for i = 1:n_groups
%     mask_order{i} = ['/Volumes/Data/' piece '/Visual/masks/' masking{i} '_mask_files.mat'];
% end