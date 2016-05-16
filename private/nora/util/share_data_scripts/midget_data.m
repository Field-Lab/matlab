% NB 2014-10-9
% get data into a reasonable shape to share

% Data Type
clear; close all;  clc
GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';%
GLMType.map_type = 'mapPRJ';
i_exp = 1; i_cell = 1;
cellselectiontype = 'all';
GLMType.fit_type = 'WN';

piece = '2012-08-09-3';
piece_file = '201208093';
stim_type = 'NSEM';
tstim = .00832750;

%% load basic info
expnumber = i_exp;
[exp_nm,cells,expname]  = cell_list( expnumber, cellselectiontype);
[StimulusPars DirPars datarun_slv datarun_mas] = Directories_Params_v23(exp_nm, GLMType.fit_type, GLMType.map_type);
clear boolean_debug map_type fit_type shead_cellID expname
inputs.exp_nm    = exp_nm;
inputs.map_type  = GLMType.map_type;
inputs.stim_type = GLMType.fit_type;

inputs.exp_nm       = exp_nm;
inputs.map_type     = GLMType.map_type;
DirPars.WN_STAdir   = NSEM_secondaryDirectories('WN_STA', inputs);
inputs.stim_type    = GLMType.fit_type;

DirPars.organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', inputs);
clear inputs

%% loading and organizing stimulus
load(['/Volumes/Lab/Users/Nora/ShareData/Data/CarlosData/' stim_type '-' piece '-StimData.mat'])
eval(['blockedmoviecell = ' stim_type 'StimData.FitMovie;'])
eval(['testmovie = ' stim_type 'StimData.TestMovie;'])
eval(['cells = fieldnames(' stim_type 'CellData);'])
clear WNStimData
% concatenate the fit movie (different blocks
n_blocks = length(blockedmoviecell);
height       = size(blockedmoviecell{1},2);
width        = size(blockedmoviecell{1},1);
fitframes    = size(blockedmoviecell{1},3);
totalframes       = n_blocks * ( fitframes) ;
concat_fullfitMovie = uint8(zeros(width, height, totalframes)) ;
for i_blk = 1:n_blocks
        framenums = ( (i_blk -1)*fitframes + 1 ) :  (i_blk *fitframes);  
        concat_fullfitMovie(:,:,framenums) = blockedmoviecell{i_blk};    
end
fitmovie = concat_fullfitMovie;
clear concat_fullfitMovie blockedmoviecell framenums height i_blk totalframes width


%%
for i_cell = 1:length(cells)
    cid = cells{i_cell};
    [~ , cell_savename, ~]  = findcelltype(cid, datarun_mas.cell_types);
    
    master_idx         = find(datarun_mas.cell_ids == cid);
    stafit_centercoord = ( datarun_mas.vision.sta_fits{master_idx}.mean );
    stafit_sd          = ( datarun_mas.vision.sta_fits{master_idx}.sd   );
    slvdim.height      = StimulusPars.slv.height; slvdim.width = StimulusPars.slv.width;
    [center_coord,sd]  = visionSTA_to_xymviCoord(stafit_centercoord, stafit_sd, StimulusPars.master, slvdim);
    clear master_idx stafit_centercoord slvdim sd
    
    eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', DirPars.organizedspikesdir, cell_savename));
    eval(sprintf('load %s/STAandROI_%s.mat STAandROI', DirPars.WN_STAdir, cell_savename));
    spikes = organizedspikes.block.t_sp_withinblock;
    STA = STAandROI.STA;
    
    % concatenate spikes
    testspikes = spikes(1:2:end);
    blockedspikes = spikes(2:2:end);
    t_start   = 0;
    T_SP = []; blk_count = 0;
    dur = tstim * fitframes;
    for k = 1:n_blocks
        blk_count = blk_count + 1;
        t_sp_full = blockedspikes{k} ; % unit of time: sec, 0 for the onset of the block
        t_sp      = t_sp_full(find(t_sp_full >  t_start));
        t_sp = t_sp - t_start;
        t_spcontext = t_sp + ( blk_count -1 )*dur;
        T_SP = [T_SP ; t_spcontext];
    end
    fitspikes = T_SP;
    clear T_SP blk_count blockedspikes k spikes t_sp t_sp_full t_spcontext t_start
    
    [~,center] = STA_Test(fitspikes, fitmovie, 1, tstim);
    fittedGLM = glm_fit(fitspikes, fitmovie, center, 'WN_STA', STA, 'monitor_refresh', 1/tstim);
    fittedGLM.xval = glm_predict(fittedGLM,testmovie, 'testspikes', testspikes);
    temp = corrcoef(conv(sum(fittedGLM.xval.rasters.glm_sim), gausswin(100)),conv(sum(fittedGLM.xval.rasters.recorded), gausswin(100)));
    fittedGLM.xval.corr = temp(2,1);
    close all
    plotfilters(fittedGLM)
    exportfig(gcf, ['/Volumes/Lab/Users/Nora/GLMFits/' piece_file '/' stim_type '/' names{cell} '_filters.eps'], 'Bounds', 'loose', 'Color', 'rgb');
    close all
    plotrasters(fittedGLM.xval, fittedGLM)
    exportfig(gcf, ['/Volumes/Lab/Users/Nora/GLMFits/' piece_file '/' stim_type '/' names{cell} '_rasters.eps'], 'Bounds', 'loose', 'Color', 'rgb');
    close all
    save(['/Volumes/Lab/Users/Nora/GLMFits/' piece_file '/' stim_type '/' names{cell} '.mat'], 'fittedGLM');

end
