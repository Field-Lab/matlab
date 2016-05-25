%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stuff to change for each piece

piece = '2016-02-17-8';

class = 'data000';
masking{1} = 'data002';
masking{2} = 'data004';
split = 'data001';
reg{1} = 'data003';
reg{2} = 'data005';

%%
% check ids across mvision dataruns

%%


% Maskin
cells{1} = [917, 2731, 3452, 4876, 5973, 7250];
n_masks{1} = 7;
sigmas{1} = [2 2 4 4 5 5 6 6 4 4 5 5 6 6];
% subgroup 1
cell_idx{1}{1} = 1:2:5;
mask_conditions{1}{1} = [1 3 5 7];
comp_conditions{1}{1} = [2 4 6 8];
% subgroup 1
cell_idx{1}{2} = 2:2:6;
mask_conditions{1}{2} = [1 9 11 13];
comp_conditions{1}{2} = [2 10 12 14];

% Amacrines
cells{2} = [2944, 4863];
n_masks{2} = 8;
sigmas{2} = [10 10 2 2 4 4 6 6 10 10 2 2 4 4 6 6];
% subgroup 1
cell_idx{2}{1} = 2;
mask_conditions{2}{1} = [1 3 5 7];
comp_conditions{2}{1} = [2 4 6 8];
% subgroup 1
cell_idx{2}{2} = 1;
mask_conditions{2}{2} = [9 11 13 15];
comp_conditions{2}{2} = [10 12 14 16];

% for i = 1:n_groups
%     mask_order{i} = ['/Volumes/Data/' piece '/Visual/masks/' masking{i} '_mask_files.mat'];
% end