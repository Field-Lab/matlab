%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stuff to change for each piece

piece = '2016-02-17-8';
Analysis_Path = ['/Volumes/Analysis/' piece '/mVision/'];

class = 'data000';
masking = 'data004';
split = 'data001';
reg{1} = 'data003';
reg{2} = 'data005';

% Maskin
cells{1} = [917, 2731, 3452, 4876, 5973, 7250];
cells{1}{1} = [917, 2731, 3452, 4876, 5973, 7250]; % these will have the same ids
cells{1}{2} = cells{1}{1};
cells{2}{1} = [918;2731;3456;4921;6065;7128]; % cell ids for concatenated mapping, matched to streamed
cells{2}{2} = cells{2}{1};
cells{3}{1} = [916;2732;3453;4923;6062;7141];% cell ids for mVision data003, matched to streamed
cells{3}{2} = [916;2733;3453;4876;5957;7144];% cell ids for mVision data005, matched to streamed
cells{4}{1} = [917;2731;3454;4877;6065;7128];% cell ids for Vision data003, matched to streamed
cells{4}{2} = [917;2733;3453;4923;5958;7129];% cell ids for Vision data005, matched to streamed

% masking params
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