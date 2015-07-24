function [StimulusPars exp_info] = Directories_Params_v23(string_piece, slv_type, map_type)
% Version 22  -- AKHeitman 2014-04-02
%             -- done 2014-04-06
% Calls:  BD = NSEM_BaseDirectories
        % NSEM_secondaryDirectories
        %StimulusPars.master = params_classificationfile(exp_info.masfile);
        %StimulusPars.slv    = params_fitandraster(slv_file);


% MFILE     Directories_Params_func  ( so that we input date and do less
%                                                      switching)
%
%
% inputs: string_date... piece number in the char form '2012-09-27-3'
%         string_fitttyepe  .. either char 'BW'   or char 'NSEM'
%         boolean_debug   either false or true,  false = long version
%                             ture is the short debugger version
%         maptype:   either  "mapEI"  or "mapPRJ"
%
%

% slv_type is an optional argument
%
% usage:  single m-file with all directories and params for the GLM
%
%
% outputs:     datarun  -loads the master classifcation (1)
%                       -loads the slave classification and spikes (2)
%              all params used throughout the GLM  
%
% AK HEITMAN START 2012-09-17 .. COMMENTING 2012-10-1  


%%%%%%%%   NOTES      %%%%%%%%%%%
% DEBUG JUST LIMITS THE NUMBER OF BLOCKS
% NSEM AND BW HAVE DIFFERENT NUMBERING SCHEME IN MATLAB STORAGE
% NOVEL - ALL NOVEL BLOCKS
% TEST  - THE STATIC RUNS
% FIT   - SUBSET OF NOVEL USED FOR MODEL FITTING
% BW AND NSEM PARAMS BASED UPON MARTIN'S STIMULUS FILES
%%
%%%%% 0.   BASE DIRECTORY THAT NEEDS TO BE SET %%%%%%

[exp_info] = experimentinfoNSEM(string_piece) ;


%% find stimulus parameters
if strcmp(slv_type , 'WN') || strcmp(slv_type, 'BW')
    slv_file = exp_info.WNfile;
elseif strcmp(slv_type, 'NSEM')
    slv_file = exp_info.NSEMfile;
end
StimulusPars.master = subR_stimparams_classificationfile(exp_info.masfile);
StimulusPars.slv    = subR_stimparams_fitandraster(slv_file);
clear slv_file
%% experiment dependent fit/test parameters

StimulusPars.slv.FitBlocks       = 6:2:(StimulusPars.slv.n_blk-2);
StimulusPars.slv.TestBlocks      = 5:2:(StimulusPars.slv.n_blk-2);
if strcmp(string_piece, '2013-10-10-0') && strcmp(slv_type, 'NSEM')
    StimulusPars.slv.n_rep = 27;
    StimulusPars.slv.n_blk = 54;
    StimulusPars.slv.NovelBlocks    =  2:2:StimulusPars.slv.n_blk; 
    StimulusPars.slv.StaticBlocks   =  1:2: (StimulusPars.slv.n_blk - 1);
    StimulusPars.slv.FitBlocks       = 4:2:54;
    StimulusPars.slv.TestBlocks      = 3:2:53;
end
StimulusPars.slv.computedtstim = .0083275;
%StimulusPars.slv.computedtstim = 1/120;



end

function params = subR_stimparams_fitandraster(string_moviename)


if strcmp(string_moviename, 'BW-8-1-0.48-11111_RNG_16807') || strcmp(string_moviename, 'WN-8-1-0.48-11111_RNG_16807')
    % Stimulus spatial temporal resolution
    params.height = 40;
    params.width =  80;
    params.frames_pertrigger = 100;
    params.tstim = .00832750;
    % params.tstim = 1/120;
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
    % params.tstim = 1/120;
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
    params.fittest_skipframes  = 120;
    params.fittest_skipseconds = 1;
    
    params.test_skipENDseconds = 60;
    params.test_skipENDframes  = 7200;
end

params.fitframes =  [params.fittest_skipframes+1 : params.frames_pernovelblock];
params.testframes = [params.fittest_skipframes+1 : params.frames_perstaticblock];


if isfield(params, 'test_skipENDframes')
    params.testframes = [(params.fittest_skipframes+1): (params.frames_perstaticblock-params.test_skipENDframes)];
end
    
% DO NOT ADD FITSECONDS OR TEST SECONDS
% THIS WON"T BE THE EXAcT NUMBER!!!
% SEONDS ARE SLIGHTLY DIFFERENT  (1/120) FROM PRECISE TIME OF SPIKES


end
function params = subR_stimparams_classificationfile(string_noisefile)


if strcmp(string_noisefile, 'RGB-10-2-0.48-11111-64x32')
    params.pixelsize = 10;
    params.height = 32; params.width  = 64;
    params.refreshrate = 2;
    params.frames_pertrigger = 50;
    params.tstim = 2 * .0083275;
    params.type = 'RGB';
    params.RNG  = 11111;
end

if strcmp(string_noisefile, 'RGB-8-1-0.48-11111-80x40')
    params.pixelsize = 8;
    params.height = 40; params.width  = 80;
    params.refreshrate = 1;
    params.frames_pertrigger = 100;
    params.tstim             = .0083275;
    params.type = 'RGB';
    params.RNG  = 11111;
    
end


end