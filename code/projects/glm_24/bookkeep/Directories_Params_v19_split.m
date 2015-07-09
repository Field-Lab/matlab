% v19_split got rid of the silly boolean_debug option   2014-01-21

% IMplicit assumption of 120 hz stimulus... perhaps need to make this more
% flexible.. we'll see.
% DirPars.NSEMprojects_home is most important!!!
% Raster:Static :: Fit:Novel

% his willcontinue to be an ugly piece of code 
% Split means we create 2 seperate dataruns .. a master and slv.
% Makes easier to see which cells faied to get mapped.
% Started 2013-12-03

%datarun_slv.block.t_frame{j} = t_frame_interpAH(datarun_slv.block.trig{j});
% Needs to be melded with GLM_AH_18
% some differences.. especially in how variables are stored
% (StimulusPars.master)


% Updatae part 0 as necessary ..everything else should stay fixed
% This version will take into account master to slv mapping types
% Different directory oranization

% get rid ofDirectories_Params all or funcall to discard any chnce of
% shadowing!   2013-11-30

% for glm_AH_18 version
% This one will introduce cleaner nomenclature for directories etc.  
% startin 2013-11-29

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

function [StimulusPars, DirPars, datarun_slv, datarun_mas] = Directories_Params_v19_split(string_date, slv_type, map_type)

%string_date = '2012-08-09-3'; boolean_debug = false; slv_type = 'BW'; nargout = 3
%if ~exist('map_type', 'var'),  map_type = 'mapEI';  end
%%%%% 0.   BASE DIRECTORY THAT NEEDS TO BE SET %%%%%%
% BD stands for base directory
BD.NSEMhome        = '/netapp/snle/lab/Experiments/Array/Analysis/akheitman/NSEM_Projects';
BD.codehome        = '/netapp/snle/lab/temp_matlabcodeAH/glm_AH_18';
% BD.baseanalysisdir = '/netapp/snle/lab/Experiments/Array/Analysis';
BD.baseanalysisdir = '/Users/colleen/Desktop/NSEM_blocked_spikes';
exp_nm = string_date;
%%%%% 0A.  EXPERIMENT , MASTER AND SLAVE, DATA INPUT  %%%
if strcmp(string_date, '2013-10-10-0'); dn_mas = 'data000'; dn_BWmapEI = 'data005';   dn_NSEMmapEI = 'data001'; end
if strcmp(string_date, '2013-08-19-6'), dn_mas = sprintf('data000'); dn_BWmapEI = sprintf('data003'); dn_NSEMmapEI = sprintf('data001'); end
if strcmp(string_date,'2012-09-27-3'), dn_mas = sprintf('data003'); dn_BWmapEI = sprintf('data005'); dn_NSEMmapEI = sprintf('data002');  end
if strcmp(string_date,'2012-08-09-3'), dn_mas = sprintf('data002'); dn_BWmapEI = sprintf('data006'); dn_NSEMmapEI = sprintf('data005');  end
if strcmp(string_date,'2012-09-21-1'), dn_mas = sprintf('data004'); dn_BWmapEI = sprintf('data006'); dn_NSEMmapEI = sprintf('data005');  end
if strcmp(string_date,'2012-09-24-2'), dn_mas = sprintf('data002'); dn_BWmapEI = sprintf('data005'); dn_NSEMmapEI = sprintf('data001');  end
if strcmp(string_date,'2012-09-13-1'), dn_mas = sprintf('data004'); dn_BWmapEI = sprintf('data007'); dn_NSEMmapEI = sprintf('data003');  end

if strcmp(string_date,'2012-04-13-4'), dn_mas = sprintf('data000'); dn_BWmapEI = sprintf('data003'); dn_NSEMmapEI = sprintf('data002');  end

if strcmp(string_date,'2012-01-27-4'), dn_mas = sprintf('data007'); dn_BWmapEI = sprintf('data003'); dn_NSEMmapEI = sprintf('data005');  end
dn_BWmapPRJ   = sprintf('%s-from-%s',   dn_BWmapEI , dn_mas );
dn_NSEMmapPRJ = sprintf('%s-from-%s', dn_NSEMmapEI , dn_mas );


%%%%% 0B.  STIMULUS PARAMETERS  %%%

%%%%% MASTER
StimulusPars.master.pixelsize = 10;
StimulusPars.master.height = 32; StimulusPars.master.width  = 64;
StimulusPars.master.refreshrate = 2;
StimulusPars.master.type = 'RGB';
if strcmp(exp_nm , '2012-08-21-1') || strcmp(exp_nm , '2012-08-09-3') || strcmp(exp_nm, '2012-04-13-4')
    StimulusPars.master.height = 40; StimulusPars.master.width  = 80; StimulusPars.master.pixelsize = 8;
end
if strcmp(exp_nm, '2012-08-09-3')
   StimulusPars.master.refreshrate = 1;
end

%%%%% BW SLV
StimulusPars.BW.commonmovie         = true;
if strcmp(exp_nm, '9999-99-99-9')
    StimulusPars.BW.commonmovie = false;
end
StimulusPars.BW.default_tstim = .00832750;
ps = 8; %pixelsize 
height = 40;
width = 80;
rr = 1; %refresh rate
fpt = 100; % frame per trigger
seedS = 11111; % stat seed 
seedA = 16807; %random number generating component A
seedM = 2^31 -1; % random number generating component M
seedC = 0  ; % rando
fit_frames = 3600;
raster_frames = 1200;

%%%%%%%%%%%%%%%
%%% NEED TO VERIFY THIS MAYBE JUST AVOID %%%
%%%%%%%%%%%%%%%
if strcmp(exp_nm, '2012-01-27-4')
    seedA = 1664525;
    seedC = 1013904223;
    seedM = 2^32;
    seedS = 11111; % initial seed
    blk = 2*[11,16,17,18,21,22,23,24,25,26,27,30,31,34,35,36,37,38,41,43,45,46,47,48,51,53,54,55,57,58];
end

if ~StimulusPars.BW.commonmovie
   display('notcommonmovieneedtheBWparams') 
end
StimulusPars.BW.height = height;
StimulusPars.BW.width = width;
StimulusPars.BW.seedRNG = 'new_seedS = mod( (seed.A*old_seedS + seedC), seedM)';
StimulusPars.BW.seedA = seedA  ; StimulusPars.BW.seedC = seedC ;
StimulusPars.BW.seedM = seedM  ; StimulusPars.BW.seedS = seedS ;
StimulusPars.BW.frames_pertrigger = fpt;
StimulusPars.BW.pixelsize = ps; StimulusPars.BW.refreshrate = rr;
StimulusPars.BW.rast_shortname = sprintf('BW-%d-%d-0.48-%d.xml' , ps , rr, seedS);
StimulusPars.BW.fit_rast_fullname = sprintf('BW-%d-%d-0.48-%d_RNG_%d', ps, rr,seedS, seedA) ;
StimulusPars.BW.ntb_o    = 13;  % TRIGGERS PER STATIC RASTER BLOCK
StimulusPars.BW.ntb_e    = 37;  % TRIGGERS PER NOVEL FITTING BLOCK    
StimulusPars.BW.nsec_o   = 10; % SECONDS PER STATIC RASTER BLOCK  
StimulusPars.BW.nsec_e   = 30; % SECONDS PER NOVEL FITTING BLOCK
StimulusPars.BW.fit_frames = fit_frames ;
StimulusPars.BW.raster_frames = raster_frames ;
clear seedA seedM seedC seedS ps rr fpt

%%%% NSEM Slv
StimulusPars.NSEM.commonmovie = true;
if strcmp(exp_nm, '9999-99-99-9')
    StimulusPars.NSEM.commonmovie = false;
end
if ~StimulusPars.NSEM.commonmovie
   display('notcommonmovieneedtheotherNSEMparams') 
end
StimulusPars.NSEM.ntb_o   = 37;  % TRIGGERS PER STATIC RASTER BLOCK
StimulusPars.NSEM.ntb_e   = 73;  % TRIGGERS PER NOVEL FITTING BLOCK    
StimulusPars.NSEM.nsec_o  = 30; % SECONDS PER STATIC RASTER BLOCK  
StimulusPars.NSEM.nsec_e  = 60; % SECONDS PER NOVEL FITTING BLOCK
StimulusPars.NSEM.default_tstim = .00832750;


%%
%%%%%%%%   1A. Fill in Dir Pars     %%%%%%%%%%
analysisdir          = sprintf('%s/%s' ,BD.baseanalysisdir,exp_nm);
DirPars.exp_nm       = string_date;
DirPars.analysisdir  = analysisdir;
DirPars.dn_mas       = dn_mas;

fitandmap    = sprintf('%s_%s', slv_type, map_type);
switch fitandmap
    case 'BW_mapEI'
        DirPars.dn_slv = dn_BWmapEI;
    case 'NSEM_mapEI'
        DirPars.dn_slv = dn_NSEMmapEI;
        case 'BW_mapPRJ'
        DirPars.dn_slv = dn_BWmapPRJ;
    case 'NSEM_mapPRJ'
        DirPars.dn_slv = dn_NSEMmapPRJ;
end

DirPars.fitandmap         = fitandmap;
DirPars.dn_BWmapEI        = dn_BWmapEI;
DirPars.dn_NSEMmapEI      = dn_NSEMmapEI;
DirPars.dn_BWmapPRJ       = dn_BWmapPRJ;
DirPars.dn_NSEMmapPRJ     = dn_NSEMmapPRJ;

DirPars.output_dir = analysisdir;

% DirPars.output_dir = sprintf('%s/GLM/%s/%s',BD.NSEMhome,exp_nm,fitandmap);
DirPars.analysisdir         = analysisdir;
DirPars.NSEMprojects_home   = BD.NSEMhome;
DirPars.stimulimatxmlfiles  = sprintf('%s/Stimuli',BD.NSEMhome);
DirPars.codedir             = BD.codehome;
DirPars.rawMovienotes       = 'NSEM rawMovie files are huge, location will vary in move from Salk to Stanford';
if ~exist(DirPars.output_dir,'dir'), mkdir(DirPars.output_dir), end
clear dn_BWmapEI dn_BWmapPRJ dn_mas dn_NSEMmapEI dn_NSEMmapPRJ

%%%% 1B. Fill in StimulusPars %%%%
StimulusPars.BW.n_rep = 60; 
StimulusPars.BW.n_blk = StimulusPars.BW.n_rep*2; 
StimulusPars.BW.NovelBlocks    =  2:2:StimulusPars.BW.n_blk; StimulusPars.BW.StaticBlocks = 1:2: (StimulusPars.BW.n_blk - 1);
StimulusPars.BW.FitBlocks      = 18:2:(StimulusPars.BW.n_blk-4); % ROUGHLY SKIP THE FIRST  MINUTES

if strcmp(exp_nm , '2012-01-27-4')
    StimulusPars.BW.FitBlocks   = 2*[11,16,17,18,21,22,23,24,25,26,27,30,31,34,35,36,37,38,41,43,45,46,47,48,51,53,54,55,57,58];
    StimulusPars.BW.NovelBlocks = 2*[11,16,17,18,21,22,23,24,25,26,27,30,31,34,35,36,37,38,41,43,45,46,47,48,51,53,54,55,57,58];
end
StimulusPars.BW.TestFrameIdx   = 1;          StimulusPars.BW.NovelFrameIdx  = 2;
StimulusPars.BW.PadVal         = 0;          

StimulusPars.BW.fr_sec = 120 ;

%%% NSEM AND BW HAVE DIFFERENT NUMBERING SCHEME IN STORAGE %%%
%%% BW FOLLOWS TRUE "STATIC THEN NOVEL" ODD EVEN REPEAT %%%
StimulusPars.NSEM.moviename = 'eye-120-3_0-3600';
StimulusPars.NSEM.frames_pertrigger = 100;
StimulusPars.NSEM.n_rep = 59;
StimulusPars.NSEM.n_blk = StimulusPars.NSEM.n_rep*2; 
StimulusPars.NSEM.fr_sec = 120 ;   

StimulusPars.NSEM.n_blk = StimulusPars.NSEM.n_rep * 2;
StimulusPars.NSEM.NovelBlocks  = 2:2:(StimulusPars.NSEM.n_blk); 
StimulusPars.NSEM.FitBlocks    = 10:2:(StimulusPars.NSEM.n_blk);
StimulusPars.NSEM.RasterBlocks = 9:2:StimulusPars.NSEM.n_blk;
StimulusPars.NSEM.StaticBlocks = 1:2: (StimulusPars.NSEM.n_blk - 1);
StimulusPars.NSEM.StaticMovieFile  = 1;
StimulusPars.NSEM.FitMovieFiles = sort( union(StimulusPars.NSEM.FitBlocks, (StimulusPars.NSEM.FitBlocks +1)), 'ascend')
StimulusPars.NSEM.NovelMovieFiles = sort( union(StimulusPars.NSEM.NovelBlocks, (StimulusPars.NSEM.NovelBlocks +1)), 'ascend')
StimulusPars.NSEM.PadVal         =.23;
%%% NSEM MOVIE PARAMS (HACK. SHOULD PULL FROM MOVIE FILE) %%%
StimulusPars.NSEM.Refresh = 8.3275;  StimulusPars.NSEM.Refreshinterp = 'millisec dur of eachframe';
StimulusPars.NSEM.width  = 80;
StimulusPars.NSEM.height = 40;

StimulusPars.NSEM.fit_frames = 7200;
StimulusPars.NSEM.raster_frames = 3600;


if strcmp(exp_nm, '2013-08-19-6')
    StimulusPars.NSEM.moviename = 'eye-long-v2';
end


% hack non-sense to see if we can get this stim to work!
if strcmp(exp_nm , '2013-10-10-0')
    StimulusPars.NSEM.moviename = 'FEM900FF_longrast';
    StimulusPars.NSEM.NovelBlocks = 2:2:54;
    StimulusPars.NSEM.FitBlocks = 2:2:54;
    StimulusPars.NSEM.RasterBlocks = 1:2:53;
    StimulusPars.NSEM.StaticBlocks = 1:2:53;
    StimulusPars.NSEM.n_rep = 27;
    StimulusPars.NSEM.n_blk = 54;
    StimulusPars.NSEM.FitMovieFiles = sort( union(StimulusPars.NSEM.FitBlocks, (StimulusPars.NSEM.FitBlocks +1)), 'ascend');
    StimulusPars.NSEM.NovelMovieFiles = sort( union(StimulusPars.NSEM.NovelBlocks, (StimulusPars.NSEM.NovelBlocks +1)), 'ascend');
    StimulusPars.NSEM.ntb_o   = 145;  % TRIGGERS PER STATIC RASTER BLOCK
    StimulusPars.NSEM.ntb_e   = 145;  % TRIGGERS PER NOVEL FITTING BLOCK    
    StimulusPars.NSEM.nsec_o  = 120; % SECONDS PER STATIC RASTER BLOCK  
    StimulusPars.NSEM.nsec_e  = 120; % SECONDS PER NOVEL FITTING BLOCK
    StimulusPars.NSEM.fit_frames = 14400;
    StimulusPars.NSEM.raster_frames = 14400;
end
    

StimulusPars.BW.evalmodel_Blocks  = StimulusPars.BW.FitBlocks - 1;  
StimulusPars.BW.RasterBlocks   = StimulusPars.BW.StaticBlocks;
StimulusPars.BW.RasterBlocks_note = 'for now this is same as static blocks but do not want to deleter yet if it is called';


StimulusPars.NSEM.evalmodel_Blocks  = StimulusPars.NSEM.FitBlocks - 1;  
StimulusPars.NSEM.RasterBlocks   = StimulusPars.NSEM.StaticBlocks;
StimulusPars.NSEM.RasterBlocks_note = 'for now this is same as static blocks but do not want to deleter yet if it is called';



%% Optional Datarun output
%%%%%%%%   2. LOAD THE DATARUN    %%%%%%%%%%%
if nargout  >= 3
	clear datarun
	if strcmp(slv_type, 'BW')
        ntb_o  =  StimulusPars.BW.ntb_o; ntb_e = StimulusPars.BW.ntb_e; 
        n_blk = StimulusPars.BW.n_blk; n_rep = StimulusPars.BW.n_rep;
	elseif strcmp(slv_type, 'NSEM')
        ntb_o  =  StimulusPars.NSEM.ntb_o; ntb_e = StimulusPars.NSEM.ntb_e;
        n_blk  = StimulusPars.NSEM.n_blk;  n_rep = StimulusPars.NSEM.n_rep;
    end
	ntb_oe = ntb_o +ntb_e;
     
    %%% SETUP DATARUN FOR LOADING THE MASTER DATA]	
%     datarun{1}.names.rrs_params_path  =[string_date, '/data003/data003.params'];
% 	datarun{1}.names.rrs_sta_path     = [string_date, '/data003/data003.sta'];

	datarun{1}.names.rrs_params_path  = sprintf('/Volumes/Analysis/%s/%s/%s.params',exp_nm,DirPars.dn_mas,DirPars.dn_mas);
	datarun{1}.names.rrs_sta_path     = sprintf('/Volumes/Analysis/%s/%s/%s.sta', exp_nm,DirPars.dn_mas,DirPars.dn_mas);
	datarun{1}.default_sta_fits       = 'vision';

	%%% SETUP DATARUN FOR ENSLAVED (BW OR NSEM .. CAN'T CLASSIFY ON OWN DUE TO MOVIE ISSUES)
%     	datarun{2}.names.rrs_neurons_path = [string_date, '/data001-from-data000/data001-from-data000.neurons'];

	datarun{2}.names.rrs_neurons_path = sprintf('/Volumes/Analysis/%s/%s-cr/%s-cr.neurons', exp_nm,DirPars.dn_slv,DirPars.dn_slv);
	if strcmp(map_type, 'mapEI')
        datarun{2}.names.map_path         = sprintf('%s/%s/cellmatch_mapEI_from_%s.txt',analysisdir,DirPars.dn_slv,DirPars.dn_mas);
    end
	datarun{2}.default_sta_fits = 'vision';

	%%% ACTUALLY LOAD UP DATA
	opt = struct('verbose',1,'load_params',1,'load_neurons',1);
	datarun = load_data(datarun,opt);
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
%         datarun_mas{1}.names.rrs_params_path  =[string_date, '/data000/data000.params'];
% 	datarun_mas{1}.names.rrs_sta_path     = [string_date, '/data000/data000.sta'];

	datarun_mas{1}.names.rrs_params_path  = sprintf('/Volumes/Analysis/%s/%s/%s.params', exp_nm, DirPars.dn_mas,DirPars.dn_mas);
	datarun_mas{1}.names.rrs_sta_path     = sprintf('/Volumes/Analysis/%s/%s/%s.sta',   exp_nm,DirPars.dn_mas,DirPars.dn_mas);
	datarun_mas{1}.default_sta_fits       = 'vision';

	%%% SETUP DATARUN FOR ENSLAVED (BW OR NSEM .. CAN'T CLASSIFY ON OWN DUE TO MOVIE ISSUES)
%     	datarun_mas{2}.names.rrs_neurons_path = [string_date, '/data001-from-data000/data001-from-data000.neurons'];
        
%     datarun_mas{2}.names.rrs_params_path  = sprintf('/Volumes/Analysis/%s/%s/%s.params', exp_nm,DirPars.dn_mas,DirPars.dn_mas);
% 	datarun_mas.names.rrs_sta_path     = sprintf('%s/%s/%s.sta',    DirPars.analysisdir,DirPars.dn_mas,DirPars.dn_mas);
	datarun_mas{2}.names.rrs_neurons_path = sprintf('/Volumes/Analysis/%s/%s-cr/%s-cr.neurons', exp_nm,DirPars.dn_slv,DirPars.dn_slv);
% 	datarun_mas.default_sta_fits       = 'vision';
    opt = struct('verbose',1,'load_params',1,'load_neurons',1);
	datarun_mas = load_data(datarun_mas,opt);

    
    
   
	%%% ADD BLOCK STRUCTURE TO THE ENSLAVED DATA
	datarun_slv.block.trig = cell(n_blk,1);
	%%% PUT TRIGGER INTO EACH BLOCK
	for k = 1:n_rep
        trg_oe = datarun_slv.triggers((k-1)*ntb_oe+1:k*ntb_oe);   %  this loop fills in the time of the triggers
        datarun_slv.block.trig{2*k-1} = trg_oe(1:ntb_o);
        datarun_slv.block.trig{2*k}   = trg_oe(ntb_o+1:end);
    end
	%keyboard
	datarun_slv.block.t_frame = cell(n_blk,1);

	for j = 1:n_blk
        datarun_slv.block.t_frame{j} = t_frame_interpAH(datarun_slv.block.trig{j});
    end

	TSTIM = nan(1,n_blk); c = 0;
	for n_blk = 1:n_blk
        c        = c+1;
        TSTIM(c) = median(diff(datarun_slv.block.t_frame{n_blk}));  
    end
	tstim               = mean(TSTIM);
	GLMPars.tstim       = tstim;
	if strcmp(slv_type, 'BW')
        StimulusPars.BW.tstim  = tstim;
    elseif strcmp(slv_type, 'NSEM')
        StimulusPars.NSEM.tstim = tstim;
    end
    %%%%%%%%%%%%%%
    if strcmp(slv_type, 'BW')
        datarun_slv.block.xml  = cell(n_blk,1);   %% empty cells for now
        sd.a = StimulusPars.BW.seedA ; 
        sd.c = StimulusPars.BW.seedC ;
        sd.m = StimulusPars.BW.seedM ; 
        sd.s = StimulusPars.BW.seedS ;
        for k = 1:n_blk     
            % KEEP TRACK OF XML FILE NAME, PUT INTO THE DATARUN 
            if rem(k,2) == 1   %%% odd blocks   should be the static movies / test movies 
                datarun_slv.block.xml{k} = sprintf('%s/Stimuli/%s/xml_files/raster/BW-8-1-0.48-11111.xml',DirPars.NSEMprojects_home,StimulusPars.BW.fit_rast_fullname);
            else
                sd.s = mod( (sd.a*sd.s + sd.c), sd.m);
                datarun_slv.block.xml{k} = sprintf('%s/Stimuli/%s/xml_files/novel/BW-8-1-0.48-%d.xml',DirPars.NSEMprojects_home,...
                    StimulusPars.BW.fit_rast_fullname , (uint32(sd.s)) );
            end    
        end
        StimulusPars.BW.xmlfiles = datarun_slv.block.xml; 
    end
end



%%

end

