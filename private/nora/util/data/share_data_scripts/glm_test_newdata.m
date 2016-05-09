piece = '2016-02-17-1';
dsave = '/Volumes/Lab/Users/Nora/Ella_Data/2016-02-17-1/';
load([dsave 'NSEMStimData_201602171.mat'],'NSEMStimData');
load([dsave 'NSEMCellData_201602171.mat'],'NSEMCellData');
names = fieldnames(NSEMCellData);
tstim = 1/NSEMStimData.monitor_refresh;

%%
cell = 201;
eval(['fitspikes = NSEMCellData.' names{cell} '.FitSpikes;'])
eval(['testspikes = NSEMCellData.' names{cell} '.TestSpikes;'])
eval(['STA = NSEMCellData.' names{cell} '.STA;'])
[~,center] = STA_Test(fitspikes, NSEMStimData.FitMovie, 1, tstim);
fittedGLM = glm_fit(fitspikes, NSEMStimData.FitMovie, center, 'WN_STA', squeeze(sum(STA,3)), 'monitor_refresh', NSEMStimData.monitor_refresh);
fittedGLM.xval = glm_predict(fittedGLM,NSEMStimData.testmovie, 'testspikes', testspikes);
figure(1); plotfilters(fittedGLM)
figure(2); plotrasters(fittedGLM.xval, fittedGLM)

