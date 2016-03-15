%% Load data
load('/Users/Nora/Desktop/research/Data/CarlosData/WN-2012-08-09-3-CellData.mat')
load('/Users/Nora/Desktop/research/Data/CarlosData/WN-2012-08-09-3-StimData.mat')

%% Get movie in shape
fitmovie = [];
for i=1:60
    fitmovie = cat(3, fitmovie, WNStimData.FitMovie{i});
end
size(fitmovie)

%% Get spikes in shape
fitspikes =[];
block_start = 0;
for i=1:60
    block_spikes = WNCellData.OFFPar_1292.Spikes{2*i}+block_start;
    fitspikes = [fitspikes; block_spikes];
    block_start = block_start+30;
end
[STA,center] = STA_Test(fitspikes, fitmovie, 1, 1/120);

%% Fitting
fittedGLM = glm_fit(fitspikes, fitmovie, center, 'WN_STA', STA);
plotfilters(fittedGLM)

%% Testing
xval = glm_predict(fittedGLM, WNStimData.TestMovie, 'testspikes', WNCellData.OFFPar_1292.Spikes(1:2:end));
plotrasters(xval, fittedGLM);

%% check out the basis we are using for the post spike filter
[~,GLMPars] = glm_parameters;
plot(prep_spikefilterbasisGP(GLMPars.spikefilters.ps, (1/120)/GLMPars.bins_per_frame));