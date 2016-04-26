% NB 2016-04-25
clear; close all;  clc

% choose cell and experiment 
exp_nm = '2013-08-19-6';
cell_savename = 'OFFPar_4174';
cell_no = 23;

% basic info
GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';%
GLMType.map_type = 'mapPRJ';
GLMType.fit_type = 'NSEM';
tstim = .00832750;
bins_per_frame = 10;

% load Ella's spikes
load('/Users/Nora/Downloads/2013-08-19-6_OFF_GLMTraces.mat')
ella_spikes =squeeze(actualSpikes(cell_no,3:end,121:end));

% load Alex's spikes
alex_fit_dir = '/Volumes/Lab/Users/akheitman/NSEM_Home/GLM_Output_Analysis/rk1_MU_PS_noCP_timekernelCONEMODEL_stimnonlin_log_powerraise/standardparams/PS_netinhibitory_domainconstrain_COB/postfilterNL_Logistic_2Par_fixMU/NSEM_mapPRJ/';
load([alex_fit_dir exp_nm '/' cell_savename '.mat'])

% load the spikes in the format I originally sent to Ella
load(['/Volumes/Lab/Users/Nora/data_files/Data/CarlosData/NSEM-' exp_nm '-CellData.mat'])
eval(['spikes = NSEMCellData.' cell_savename '.Spikes;'])

%% Step one: Can I make the spikes I sent to Ella look like Alex's?

% load sitmulus parameters, including how many seconds alex skips, etc
StimulusPars = Directories_Params_v23(exp_nm, GLMType.fit_type, GLMType.map_type);
rasterblocks = StimulusPars.slv.TestBlocks;
t_start      = StimulusPars.slv.fittest_skipseconds;

% sort the spikes into cells, ignoring the initial skipped time
raster_spiketimes = cell(length(rasterblocks),1);
for i_blk = 1 : length(rasterblocks)
	blknum = rasterblocks(i_blk);
	sptimes = spikes{blknum} - t_start;
	sptimes = sptimes(find(sptimes > 0 ) );
    if isfield(StimulusPars.slv, 'test_skipENDseconds')
        sptimes = sptimes(find(sptimes < (StimulusPars.slv.test_skipENDseconds - StimulusPars.slv.fittest_skipseconds - .1)));
    end
    raster_spiketimes{i_blk} = sptimes;
end 

% sort the cells into logical array 
params.bins       = bins_per_frame *length(StimulusPars.slv.testframes); 
params.trials     = length(raster_spiketimes); 
params.bindur = tstim/bins_per_frame;
logicalspike = zeros(params.trials,params.bins) ;         
for i_blk = 1 : params.trials
    spt = raster_spiketimes{i_blk};
    binnumber = ceil(spt / params.bindur );
    logicalspike( i_blk, binnumber )  =  logicalspike( i_blk,binnumber ) + 1;
end 
clear i_blk spt sptimes

% compare this logical array with alex's stored spikes, output should be 0
sum(any(fittedGLM.xvalperformance.rasters.recorded - logicalspike))

%% Okay that looks good. 
% I think the issue is the t_start = StimulusPars.slv.fittest_skipseconds = 1 on line 29 above. 
% So can we make Alex's spikes look like Ella's using 0.993 instead of 1?

% change t_start
t_start      = 120*tstim;

% sort the spikes into cells, ignoring the initial skipped time
raster_spiketimes = cell(length(rasterblocks),1);
for i_blk = 1 : length(rasterblocks)
	blknum = rasterblocks(i_blk);
	sptimes = spikes{blknum} - t_start;
	sptimes = sptimes(find(sptimes > 0 ) );
    if isfield(StimulusPars.slv, 'test_skipENDseconds')
        sptimes = sptimes(find(sptimes < (StimulusPars.slv.test_skipENDseconds - StimulusPars.slv.fittest_skipseconds - .1)));
    end
    raster_spiketimes{i_blk} = sptimes;
end 

% sort the cells into logical array
bins_per_frame = 1;
params.bins       = bins_per_frame *length(StimulusPars.slv.testframes);
params.trials     = length(raster_spiketimes);
params.bindur = tstim/bins_per_frame;
logicalspike = zeros(params.trials,params.bins) ;
for i_blk = 1 : params.trials
    spt = raster_spiketimes{i_blk};
    for i = 1:length(spt)
        binnumber = ceil(spt(i) / params.bindur );
        logicalspike( i_blk, binnumber )  =  logicalspike( i_blk,binnumber ) + 1;
    end
end

% compare this logical array with ella's stored spikes, output should be 0
sum(any(ella_spikes - logicalspike))

