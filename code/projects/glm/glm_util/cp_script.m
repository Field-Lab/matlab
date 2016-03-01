clear;

exp = '2012-08-09-3';
cell = 'ONPar_841';
fitname = 'rk2_MU_PS_CP_p8IDp8';
type = 'NSEM';

load(['/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/good_fits/' fitname '/standardparams/' type '_mapPRJ/' exp '/' cell '.mat']);
secondDir.exp_nm    = exp;
secondDir.map_type  = 'mapPRJ';
secondDir.stim_type = type;
secondDir.fitname   = fittedGLM.GLMType.fitname;
Dirs.organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', secondDir);
% loading the neighboring spikes to neighborspikes.home
for i_pair=1:fittedGLM.GLMPars.spikefilters.cp.n_couplings
    eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', Dirs.organizedspikesdir,  fittedGLM.cellinfo.pair_savename{i_pair}));
    neighborspikes{i_pair} = organizedspikes.block.t_sp_withinblock(1:2:end);
end

testmovie0 = loadmoviematfile(exp , type, fittedGLM.GLMType.cone_model,'testmovie');
xval = glm_predict(fittedGLM, testmovie0{1}.matrix, 'neighborspikes', neighborspikes);

plotrasters(fittedGLM.xvalperformance, fittedGLM);