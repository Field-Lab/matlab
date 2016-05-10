function GLM_PSTH_scatter(exp_nm, cid, GLMType)

% Load some parameters
GLMType.fit_type='NSEM';
[StimulusPars, DirPars, ~, datarun_mas] = Directories_Params_v23(exp_nm, GLMType.fit_type, GLMType.map_type);
[~ , cell_savename, ~]  = findcelltype(cid, datarun_mas.cell_types);
inputs.exp_nm    = exp_nm;
inputs.map_type  = GLMType.map_type;
GLMType.fitname  = GLM_fitname(GLMType);
inputs.fitname   = GLMType.fitname;

% Load the white noise GLM fit
GLMType.fit_type='WN';
inputs.stim_type = GLMType.fit_type;
d_save = NSEM_secondaryDirectories('savedir_GLMfit', inputs);
load([d_save '/' cell_savename '.mat'])
fittedGLM_WN=fittedGLM;
fittedGLM_WN.GLMType.Subunits=false;

% Load the NSEM GLM fit
GLMType.fit_type = 'NSEM';
inputs.stim_type = GLMType.fit_type;
d_save = NSEM_secondaryDirectories('savedir_GLMfit', inputs);
load([d_save '/' cell_savename '.mat'])
fittedGLM_NSEM=fittedGLM;

% Load the NSEM stimulus
[testmovie] = loadmoviematfile(exp_nm , GLMType.fit_type, GLMType.cone_model,'testmovie');

% Load the spikes
if GLMType.CBP
   DirPars.organizedspikesdir = ['/Volumes/Analysis/nora/NSEM/CBPBlockedSpikes/' exp_nm '/NSEM_mapPRJ'];
else
    DirPars.organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', inputs);
end
clear inputs

eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', DirPars.organizedspikesdir, cell_savename));
cell_organizedspikes=organizedspikes;

if GLMType.CouplingFilters
    pairs=fittedGLM_NSEM.cellinfo.pairs;
    n_couplings=length(pairs);
    for j=1:n_couplings
        [~ , pair_savename{j}, ~]  = findcelltype(pairs(j), datarun_mas.cell_types);
        eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', DirPars.organizedspikesdir, pair_savename{j}));
        neighbor_organizedspikes{j}=organizedspikes;
    end
else
    neighbor_organizedspikes=0;
end

% Performance and plotting
xvalperformance_cross = eval_xvalperformance_NEW_CP(fittedGLM_WN, StimulusPars.slv, cell_organizedspikes,neighbor_organizedspikes,testmovie);
plotscatter(xvalperformance_cross,fittedGLM_NSEM)
plotscatter(fittedGLM_NSEM.xvalperformance,fittedGLM_NSEM)
plotscatter(fittedGLM_WN.xvalperformance,fittedGLM_WN)

end