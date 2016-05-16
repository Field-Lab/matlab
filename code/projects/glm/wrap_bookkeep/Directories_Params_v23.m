function [StimulusPars, DirPars, datarun_slv, datarun_mas] = Directories_Params_v23(string_piece, slv_type, map_type)
% JAVA CODE %%
%%% ACTUALLY LOAD UP DATA   JAVA CALLS
%	opt = struct('verbose',1,'load_params',1,'load_neurons',1);%
%	datarun = load_data(datarun,opt);
%	if isfield(datarun{2}.names,'map_path')
%        datarun=load_map(datarun);
%    else
%        datarun=map_cell_types(datarun,'verbose',true);
%    end

% subR stands for subRoutine

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

[BD] = NSEM_BaseDirectories;
if ~exist('map_type', 'var'),  map_type = 'mapPRJ';  end
[exp_info] = experimentinfoNSEM(string_piece) ;
if strcmp(slv_type , 'WN') || strcmp(slv_type, 'BW')
    dn_slv = exp_info.dr.slvWN;
elseif strcmp(slv_type, 'NSEM')
    dn_slv = exp_info.dr.slvNSEM;
end
if strcmp(map_type, 'mapPRJ')
    dn_slv  = sprintf('%s-from-%s', dn_slv , exp_info.dr.mas);
end

DirPars.exp_nm       = string_piece;
DirPars.slv_type     = slv_type;
DirPars.map_type     = map_type;
DirPars.dn_mas       = exp_info.dr.mas;
DirPars.dn_slv       = dn_slv;

DirPars.analysisdir     = sprintf('%s/%s' ,BD.analysisdir,string_piece);

%{
inputs.stim_type  = slv_type;
inputs.map_type   = map_type;
inputs.exp_nm     = string_piece;
directory_type    = 'GLM_fitandraster';
DirPars.output_dir   = NSEM_secondaryDirectories(directory_type, inputs);
%}

DirPars.codedir             = BD.GLM_codehome;
clear inputs analysisdir directory_type 


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


%% Optional Datarun output
%%%%%%%%   2. LOAD THE DATARUN    %%%%%%%%%%%
SPars = StimulusPars.slv;
if nargout  >= 3
	clear datarun
	if strcmp(slv_type, 'BW') || strcmp( slv_type, 'WN')
        ntb_static  =  SPars.triggers_perstaticblock; ntb_novel = SPars.triggers_pernovelblock; 
        n_blk = SPars.n_blk; n_rep = SPars.n_rep;
	elseif strcmp(slv_type, 'NSEM')
        ntb_static  =  SPars.triggers_perstaticblock; ntb_novel = SPars.triggers_pernovelblock; 
        n_blk  = SPars.n_blk;  n_rep = SPars.n_rep;
    end
	ntb_combined = ntb_static +ntb_novel;
     
    %%% SETUP DATARUN FOR LOADING THE MASTER DATA
    %{
	datarun{1}.names.rrs_params_path  = sprintf('%s/%s/%s.params', DirPars.analysisdir,DirPars.dn_mas,DirPars.dn_mas);
	datarun{1}.names.rrs_sta_path     = sprintf('%s/%s/%s.sta',    DirPars.analysisdir,DirPars.dn_mas,DirPars.dn_mas);
%}		
    opt = struct('verbose',1,'load_params',1,'load_neurons',1);
    datarun{1} = load_data([DirPars.analysisdir '/' DirPars.dn_mas '/' DirPars.dn_mas], opt);
    datarun{1}.default_sta_fits       = 'vision';

	%%% SETUP DATARUN FOR ENSLAVED (BW OR NSEM .. CAN'T CLASSIFY ON OWN DUE TO MOVIE ISSUES)
	%datarun{2}.names.rrs_neurons_path = sprintf('%s/%s/%s.neurons', DirPars.analysisdir,DirPars.dn_slv,DirPars.dn_slv);
	datarun{2} = load_data([DirPars.analysisdir '/' DirPars.dn_slv '/' DirPars.dn_slv], opt);
    if strcmp(map_type, 'mapEI')
        datarun{2}.names.map_path         = sprintf('%s/%s/cellmatch_mapEI_from_%s.txt',analysisdir,DirPars.dn_slv,DirPars.dn_mas);
    end
	datarun{2}.default_sta_fits = 'vision';

	%%% ACTUALLY LOAD UP DATA   JAVA CALLS
    %datarun = load_data(datarun,opt);
	if isfield(datarun{2}.names,'map_path')
        datarun=load_map(datarun);
    else
        datarun=map_cell_types(datarun,'verbose',true);
    end
    clear opt
    
    datarun_slv = datarun{2};
    clear datarun
    
    %%% RELOAD MASTER WITH FULL 
    %%% HACK TO GET AROUND WEIRD CELL DROPPING IN LOAD_MAP or
    %%% MAP_CELL_TYPES
    datarun_mas.names.rrs_params_path  = sprintf('%s/%s/%s.params', DirPars.analysisdir,DirPars.dn_mas,DirPars.dn_mas);
	datarun_mas.names.rrs_sta_path     = sprintf('%s/%s/%s.sta',    DirPars.analysisdir,DirPars.dn_mas,DirPars.dn_mas);
    datarun_mas.names.rrs_neurons_path = sprintf('%s/%s/%s.neurons', DirPars.analysisdir,DirPars.dn_mas,DirPars.dn_mas);
	datarun_mas.default_sta_fits       = 'vision';
    opt = struct('verbose',1,'load_params',1,'load_neurons',1);
	datarun_mas = load_data(datarun_mas,opt);

    
    
   
	%%% ADD BLOCK STRUCTURE TO THE ENSLAVED DATA
	datarun_slv.block.trig = cell(n_blk,1);
	%%% PUT TRIGGER INTO EACH BLOCK
	for k = 1: n_rep
        trg_oe = datarun_slv.triggers((k-1)*ntb_combined+1:k*ntb_combined);   %  this loop fills in the time of the triggers
        datarun_slv.block.trig{2*k-1} = trg_oe(1:ntb_static);
        datarun_slv.block.trig{2*k}   = trg_oe(ntb_static+1:end);
    end
	%keyboard
	datarun_slv.block.t_frame = cell(n_blk,1);

	for j = 1:n_blk
        datarun_slv.block.t_frame{j} = subR_t_frame_interpAH(datarun_slv.block.trig{j});
    end

	TSTIM = nan(1,n_blk); c = 0;
	for n_blk = 1:n_blk
        c        = c+1;
        TSTIM(c) = median(diff(datarun_slv.block.t_frame{n_blk}));  
    end
	tstim               = mean(TSTIM);
    SPars.computedtstim = tstim;
	if strcmp(slv_type, 'BW') || strcmp( slv_type, 'WN')
        SPars.tstim  = tstim;
    elseif strcmp(slv_type, 'NSEM')
        SPars.tstim = tstim;
    end
    %%%%%%%%%%%%%%
    if strcmp(slv_type, 'BW') || strcmp( slv_type, 'WN')
        datarun_slv.block.xml  = cell(n_blk,1);   %% empty cells for now
        sd.a = SPars.seedA ; 
        sd.c = SPars.seedC ;
        sd.m = SPars.seedM ; 
        sd.s = SPars.seedS ;
        for k = 1:n_blk     
            % KEEP TRACK OF XML FILE NAME, PUT INTO THE DATARUN 
            if rem(k,2) == 1   %%% odd blocks   should be the static movies / test movies 
                datarun_slv.block.xml{k} = sprintf('%s/Stimuli/%s/xml_files/raster/BW-8-1-0.48-11111.xml',BD.NSEM_home,SPars.novel_static_fullname);
            else
                sd.s = mod( (sd.a*sd.s + sd.c), sd.m);
                datarun_slv.block.xml{k} = sprintf('%s/Stimuli/%s/xml_files/novel/BW-8-1-0.48-%d.xml',BD.NSEM_home,...
                    SPars.novel_static_fullname , (uint32(sd.s)) );
            end    
        end
        SPars.xmlfiles = datarun_slv.block.xml; 
    end
    StimulusPars.slv.computedtstim = tstim;
end

end

function t_frame = subR_t_frame_interpAH(t_trig,interp_method,fl_fig)
% return the timing of stimulus frames, based on 100-frames trigger.
% Usage: t_frame = t_frame_interp(t_trig,interp_method,fl_fig)

% AKHeitman  inheriting and modifying 2013-12-04
% edoi@salk.edu, 2011-11-16.

if ~exist('interp_method','var')
    interp_method = 'linear';
end
if ~exist('fl_fig','var')
    fl_fig = 0;
end

N = 100; % number of frames per trigger
n_frame = (length(t_trig)-1)*N+1; % total no of frames
idx_frame = 0:(n_frame-1);   % subtract 1 because the index starts with 0.
idx_trig  = 0:N:(n_frame-1); % ditto.
if idx_frame(end) ~= idx_trig(end), fprintf('error\n'), end

t_frame = interp1(idx_trig,t_trig,idx_frame,interp_method); 
t_frame = t_frame(1:(end-1));

% check with visualization 
if fl_fig
    figure(1), clf
    plot(idx_trig,t_trig,'r.');
    hold on
    plot(idx_frame(1:(end-1)),t_frame,'b')
    xlabel(sprintf('frame idx [0,%d]',idx_trig(end)))
    xlim([0,idx_trig(end)])
    ylabel('time [sec]')
    title('Frame time')
    
    figure(2), clf
    plot(idx_frame(1:(end-2)),diff(t_frame),'b.')
    xlim([0,idx_trig(end-1)])
    xlabel(sprintf('frame idx [0,%d]',idx_trig(end)-1))
    ylabel('time [sec]')
    title('Difference of frame time')

    figure(3), clf
    hist(diff(t_trig))
    title('Histogram of difference of triggers')
    % the jitter is likely to be in the unit of sampling interval, i.e.,
    % 1/datarun.sampling_rate.
end
end

function params = subR_stimparams_fitandraster(string_moviename)


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