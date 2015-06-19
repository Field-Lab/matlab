datapath = '/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk1_MU_PS_noCP_SU_p8IDp8/standardparams/WN_mapPRJ/2012-08-09-3/';
matfiles = dir([datapath 'ONPar*.mat']);

stim = zeros(80);
stimsize.width = 80;
stimsize.height = 40;
length = 11;

for file = 1:size(matfiles,1)
   load([datapath matfiles(file).name]);
   ROI = ROI_coord(length, fittedGLM.cellinfo.slave_centercoord, stimsize);
   stim(ROI.yvals, ROI.xvals) = abs(fittedGLM.linearfilters.Stimulus.space_rk1);
   disp(file)
end


%%
clear datapath
datapath{1} = '/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk1_MU_PS_noCP_SU_p8IDp8/standardparams/WN_mapPRJ/2012-08-09-3/';
datapath{2} = '/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk1_MU_PS_noCP_p8IDp8/standardparams/WN_mapPRJ/2012-08-09-3/';

cid = zeros(40,2);
BPS = zeros(40,2);

for fit = 1:2
    matfiles = dir([datapath{fit} 'ONPar*.mat']);
    
    for file = 1:size(matfiles,1)
        load([datapath{fit} matfiles(file).name]);
        BPS(file,fit) = fittedGLM.xvalperformance.glm_normedbits;
        cid(file,fit) = fittedGLM.cellinfo.cid;
    end
    
end

%%
scatter(BPS(:,1), BPS(:,2))
hold on; plot([0 1], [0 1])
axis square
