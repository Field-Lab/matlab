%AKHeitman 2014-04-03

%%%% CAUTION %%%%%%%%%%%
% DO NOT ADD FITSECONDS OR TEST SECONDS
% THIS WON"T BE THE EXAcT NUMBER!!!
% SEONDS ARE SLIGHTLY DIFFERENT  (1/120) FROM PRECISE TIME OF SPIKES
%%%%%%%%%%%%%%%%%%%%%%%%%

%Strange Caveat.. 2012-01-27-4
%basically ignor this data set .. it is awful to begin with
%{
if strcmp(exp_nm, '2012-01-27-4')
    seedA = 1664525;
    seedC = 1013904223;
    seedM = 2^32;
    seedS = 11111; % initial seed
    blk = 2*[11,16,17,18,21,22,23,24,25,26,27,30,31,34,35,36,37,38,41,43,45,46,47,48,51,53,54,55,57,58];
end
%}

% this can all be tested by looking if we generate good rasters from here!

function params = stimparams_fitandraster(string_moviename)


if strcmp(string_moviename, 'BW-8-1-0.48-11111_RNG_16807') || strcmp(string_moviename, 'WN-8-1-0.48-11111_RNG_16807')
    % Stimulus spatial temporal resolution
    params.height = 40;
    params.width =  80;
    params.frames_pertrigger = 100;
    params.tstim = .00832750;
    params.pixelsize = 8; 
    params.refreshrate = 1;
    params.fr_sec = 120 ;
    params.avgIntensity = .5;

    % Fit / Raster structure 
    params.triggers_perstaticblock    = 13;  % TRIGGERS PER STATIC RASTER BLOCK
    params.triggers_pernovelblock   = 37;  % TRIGGERS PER NOVEL FITTING BLOCK    
    params.seconds_perstaticblock   = 10; % SECONDS PER STATIC RASTER BLOCK  
    params.seconds_pernovelblock  = 30; % SECONDS PER NOVEL FITTING BLOCK
    params.frames_pernovelblock = 3600;
    params.frames_perstaticblock = 1200;
    params.n_rep = 60; 
    params.n_blk = params.n_rep*2; 
    params.NovelBlocks    =  2:2:params.n_blk; 
    params.StaticBlocks   =  1:2: (params.n_blk - 1);
    
    % Fitting and testing time we should ignore due to jumps in stimulus
    params.fittest_skipframes = 60;
    params.fittest_skipseconds = .5;

    % White noise specific parameters
    params.static_shortname = sprintf('BW-8-1-0.48-11111.xml');
    params.novel_static_fullname = string_moviename;
    params.seedRNG = 'new_seedS = mod( (seed.A*old_seedS + seedC), seedM)';
    seedS = 11111; % stat seed 
    seedA = 16807; %random number generating component A
    seedM = 2^31 -1; % random number generating component M
    seedC = 0  ; % rando
    params.seedA = seedA  ; params.seedC = seedC ;
    params.seedM = seedM  ; params.seedS = seedS ;
end

if strcmp(string_moviename, 'NSEM_eye-120-3_0-3600') || strcmp(string_moviename, 'NSEM_eye-long-v2')
    % Stimulus spatial temporal resolution
    params.height = 40;
    params.width =  80;
    params.frames_pertrigger = 100;
    params.tstim = .00832750;
    params.pixelsize = 8; 
    params.refreshrate = 1;
    params.fr_sec = 120 ;
    params.avgIntensity = .23;
    
    % Fit / test structure 
    params.triggers_perstaticblock   = 37;  % TRIGGERS PER STATIC RASTER BLOCK
    params.triggers_pernovelblock  = 73;  % TRIGGERS PER NOVEL FITTING BLOCK    
    params.seconds_perstaticblock  = 30; % SECONDS PER STATIC RASTER BLOCK  
    params.seconds_pernovelblock = 60; % SECONDS PER NOVEL FITTING BLOCK
    params.frames_perstaticblock = 3600;
	params.frames_pernovelblock  = 7200;
    params.tstim = .00832750;
    params.n_rep = 60; 
    params.n_blk = params.n_rep*2; 
    
    params.NovelBlocks    =  2:2:params.n_blk; 
    params.StaticBlocks   =  1:2: (params.n_blk - 1);
    % Fitting and testing time we should ignore due to jumps in stimulus
    params.fittest_skipframes = 120;
    params.fittest_skipseconds = 1;
end

if strcmp(string_moviename, 'NSEM_FEM900FF_longrast')
    % Stimulus spatial temporal resolution
    params.height = 40;
    params.width =  80;
    params.frames_pertrigger = 100;
    params.tstim = .00832750;
    params.pixelsize = 8; 
    params.refreshrate = 1;
    params.fr_sec = 120;
    params.avgIntensity = .23;
    
    % Fit / Raster structure   
    params.triggers_perstaticblock   = 145;  % TRIGGERS PER STATIC RASTER BLOCK
    params.triggers_pernovelblock  = 145;  % TRIGGERS PER NOVEL FITTING BLOCK    
    params.seconds_perstaticblock  = 120; % SECONDS PER STATIC RASTER BLOCK  
    params.seconds_pernovelblock = 120; % SECONDS PER NOVEL FITTING BLOCK
    params.frames_pernovelblock = 14400;
    params.frames_perstaticblock = 14400;
    params.tstim = .00832750;
    params.n_rep = 30;
    params.n_blk = 60;
    params.NovelBlocks    =  2:2:params.n_blk; 
    params.StaticBlocks   =  1:2: (params.n_blk - 1);
    params.fittest_skipframes = 120;
    params.fittest_skipseconds = 1;
end

params.fitframes =  [params.fittest_skipframes+1 : params.frames_pernovelblock];
params.testframes = [params.fittest_skipframes+1 : params.frames_perstaticblock];
% DO NOT ADD FITSECONDS OR TEST SECONDS
% THIS WON"T BE THE EXAcT NUMBER!!!
% SEONDS ARE SLIGHTLY DIFFERENT  (1/120) FROM PRECISE TIME OF SPIKES


end
