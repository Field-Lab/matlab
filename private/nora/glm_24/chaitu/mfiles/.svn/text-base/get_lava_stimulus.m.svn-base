function [X stim_params contrast carrier state] = get_lava_stimulus(sz,samples)


% generate movie with correct parameters

% rf size
%sz = 4 ;

% file names
%file_obv = ['~/data/local/stim/lava-' num2str(sz) '.rawMovie']
%file_obv = ['./lava-temp.rawMovie'];
%file_avi = ['./lava-temp.avi'];

% frames per second
%fps = 120;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% make parameters

T_band = load('~/data/lava/T_band.mat');
T_con = load('~/data/lava/T_con.mat');
T_band = T_band.K_band;
T_con = T_con.y';

stim_params = struct('rf_size',sz, 'rf_surround',4, 'seed', 11111, ...
		     'dims', [80 160 samples], 'bandFilt', T_band, ...
		     'conFilt', T_con);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% generate movie

%[stim_vec, contrast_vec] = lava(stim_params, [], state);
[X contrast carrier state] = lava(stim_params, []);