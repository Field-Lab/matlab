%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stuff to change for each piece

piece = '2016-02-17-8';
Analysis_Path = ['/Volumes/Analysis/' piece '/mVision/'];
fig_save = ['/Volumes/Lab/Users/Nora/Fig_Output/' piece '/mVision/'];


masking = 'data002';
split = 'data001';
reg{1} = 'data003';
reg{2} = 'data005';
class = ['/Volumes/Analysis/' piece '/streamed/data000/data000'];

cells_orig = [917, 2731, 3452, 4876, 5973, 7250];

%{
cells_reg{1}= [917, 2731, 3452, 4876, 5973, 7250];
cells_reg{2}= [917, 2731, 3452, 4876, 5973, 7250];
cells_masking= [917, 2731, 3452, 4876, 5973, 7250];
%}

%%{
cells_reg{1} = [916;2732;3453;4923;6062;7141];% cell ids for mVision data003, matched to streamed
cells_reg{2} = [916;2733;3453;4876;5957;7144];% cell ids for mVision data005, matched to streamed
cells_masking = [918;2732;3471;4877;5962;7143];% cell ids for mVision data002, matched to streamed
%}

% masking params
n_masks = 7;
sigmas = [2 2 4 4 5 5 6 6 4 4 5 5 6 6];
% subgroup 1
cell_idx{1} = 1:2:5;
mask_conditions{1} = [1 3 5 7];
comp_conditions{1} = [2 4 6 8];
% subgroup 1
cell_idx{2} = 2:2:6;
mask_conditions{2} = [1 9 11 13];
comp_conditions{2} = [2 10 12 14];