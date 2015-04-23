
% Fill in whatever we forgot to put in each structure

clear; clc; close all
%computation_type = 'rasterprecision_paireddistance_Vector';
%computation_type = 'rasterprecision_paireddistance_Viktor';
computation_type = 'rasterprecision_prediction_avgsignal';
exps = 1:4; %i_exp = 4;
stimtypes = 1:2; %i_stimtype = 2;
celltypes = 1:2; %i_celltype = 1; i_cell =1;
%}

%%

% SET DIRECTORIES
BD = NSEM_BaseDirectories;
Dirs.rastmetdir = sprintf('%s/raster_performance', BD.Cell_Selection);
if ~exist(Dirs.rastmetdir, 'dir'), mkdir(Dirs.rastmetdir); end
Dirs.rast_dir = BD.BlockedSpikes;

% LOAD CELL NUMBERS / INITIALIZE SOLUTION STRUCTURE
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection)); 
raster_scores = allcells;

% HARD PARAMETERS OF THE RASTER WE ALREADY HAVE SET
hard_params.raster_params.bindur         = .00083275;
hard_params.raster_params.bins_per_frame = 10;
hard_params.map_type = 'mapPRJ';

% TIME SCALES AND PARAMETERS
if strcmp(computation_type, 'rasterprecision_prediction_avgsignal')
    metric_params.smoothbins = [.5 1 2 4 7 10 13 17 20 25 30 40 50 65 80 100 120 150 200 250 350 500 1000];
end
if strcmp(computation_type, 'rasterprecision_paireddistance_Viktor')
    metric_params.ViktorTimeParams_Bins = [1 2 4 8 16 32 64 128 256 512 1024];
    metric_params.pairnumbers           = 50;
    metric_params.bindur                = hard_params.raster_params.bindur;
    
    %metric_params.ViktorTimeParams_Bins    = [4 8 16 32 64 128 256];% 512 1024];
    %metric_params.pairnumbers	= 30
end
if strcmp(computation_type, 'rasterprecision_paireddistance_Vector')
   metric_params.smoothbins    = [1 2 4 8 16 32 64 128 256 512 1024];
   metric_params.pairnumbers	= 50; 
end
    
% HANDLE VARARGIN/OPTIONS
%{
argnum = nargin
if argnum == 4
    display( 'No special changes to Computation: Raw Computation on Unaltered Rasters')
end
if argnum > 4
    if ~strcmp(varargin{1}.option_type , 'Error')
    end
end
%}



%% Loop experiments/stimulus/celltypes to execute raster compuation
for i_exp = exps
    for i_stimtype = stimtypes   
        for i_celltype = celltypes
            %% Experiment Dependent Parameters
            
            % CLEAN UP
            clear StimPars raster_params secondDir
            % LOAD STIMULUS PARAMETERS / DEFINE CELL NUMBERS
            exp_nm  = allcells{i_exp}.exp_nm;
            expname = allcells{i_exp}.expname;
            map_type= hard_params.map_type;
            if i_stimtype == 1, stimtype = 'WN';   end
            if i_stimtype == 2, stimtype = 'NSEM'; end
            if i_celltype == 1; cellgroup = allcells{i_exp}.ONP;  celltype = 'ONPar'; end
            if i_celltype == 2; cellgroup = allcells{i_exp}.OFFP; celltype = 'OFFPar'; end
            [StimPars]    = Directories_Params_v23(exp_nm, stimtype, map_type);
            
            % INCORPORATE STIMULUS PARAMETERS INTO RASTER_PARAMS
            raster_params                       = hard_params.raster_params;
            raster_params.evalblocks            = StimPars.slv.TestBlocks;
            raster_params.bins                  = raster_params.bins_per_frame *length(StimPars.slv.testframes);
            raster_params.fittest_skipseconds   = StimPars.slv.fittest_skipseconds;
            
            % TAKE CARE OF DIRECTORIES
            secondDir.exp_nm = exp_nm;
            secondDir.stim_type = stimtype;
            secondDir.map_type  = 'mapPRJ';
            Dirs.organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', secondDir);
            Dirs.savedir = sprintf('%s/%s', Dirs.organizedspikesdir,computation_type);
            %if ~exist(Dirs.savedir,'dir'), mkdir(Dirs.savedir); end
            
            %cellgroup = fliplr(cellgroup)
            %%
            for i_cell = 1:length(cellgroup)
                
                % CLEAN UP
                clear raster_metrics organizedspikes raster_precision
                
                % EXTRACT BINNED RASTER
                cid = cellgroup(i_cell);
                cell_savename = sprintf('%s_%d', celltype,cid);
                 display(sprintf('Working on Exp: %s Stim: %s Cell: %s', exp_nm,stimtype, cell_savename));              
                eval(sprintf('load %s/%s_%s.mat raster_metrics', Dirs.savedir,computation_type, cell_savename));
                raster_metrics.hard_params = hard_params;
                if strcmp(computation_type, 'rasterprecision_paireddistance_Vector')
                    raster_metrics.paireddistances_vector.smoothing_bins      = metric_params.smoothbins;
                end
                eval(sprintf('save %s/%s_%s.mat raster_metrics', Dirs.savedir,computation_type, cell_savename));
            end
            
            
            
            
        end
    end
end