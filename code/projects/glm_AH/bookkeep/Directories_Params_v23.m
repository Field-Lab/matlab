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

%{
clear all
string_piece = '2012-08-09-3'; slv_type = 'NSEM'; map_type = 'mapPRJ'; 
Directories_Params_v22(string_piece, slv_type, map_type)
[StimulusPars DirPars datarun_slv datarun_mas] = Directories_Params_v22(string_piece, slv_type, map_type)
%}
function [StimulusPars DirPars datarun_slv datarun_mas] = Directories_Params_v23(string_piece, slv_type, map_type)

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
StimulusPars.master = stimparams_classificationfile(exp_info.masfile);
StimulusPars.slv    = stimparams_fitandraster(slv_file);
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
	datarun{1}.names.rrs_params_path  = sprintf('%s/%s/%s.params', DirPars.analysisdir,DirPars.dn_mas,DirPars.dn_mas);
	datarun{1}.names.rrs_sta_path     = sprintf('%s/%s/%s.sta',    DirPars.analysisdir,DirPars.dn_mas,DirPars.dn_mas);
	datarun{1}.default_sta_fits       = 'vision';

	%%% SETUP DATARUN FOR ENSLAVED (BW OR NSEM .. CAN'T CLASSIFY ON OWN DUE TO MOVIE ISSUES)
	datarun{2}.names.rrs_neurons_path = sprintf('%s/%s/%s.neurons', DirPars.analysisdir,DirPars.dn_slv,DirPars.dn_slv);
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
        datarun_slv.block.t_frame{j} = t_frame_interpAH(datarun_slv.block.trig{j});
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




%%

end

