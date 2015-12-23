clear all

neuronID = 212;
stimElec = 22;
offsets = [0 10 20 40 80 160];
pathToElecResp = '/snle/lab/Experiments/Array/Analysis/2011-08-01-0/data008/';
ampRange = [0.8 2.1]; %used for plotting thresholds that are out of range and can't be determined
threshYLims = [0.8 2.1];


for x = 1; %for code folding purposes only
    prePulse(1).elec = 4;
    prePulse(1).amps = 2.510;
    prePulse(1).analyzedPreResponse = 1; %1: analyzed, 2: not analyzed because almost all successes, 0: not analyzed for other reason
    
    prePulse(2).elec = 12;
    prePulse(2).amps = 2.510;
    prePulse(2).analyzedPreResponse = 1;
    
    prePulse(3).elec = 14;
    prePulse(3).amps = 2.510;
    prePulse(3).analyzedPreResponse = 1;
    
    prePulse(4).elec = 15;
    prePulse(4).amps = [2.209 2.409 2.510];
    prePulse(4).analyzedPreResponse = [1 1 1];
    
    prePulse(5).elec = 17;
    prePulse(5).amps = [2.209 2.409 2.510];
    prePulse(5).analyzedPreResponse = [1 1 1];
    
    prePulse(6).elec = 18;
    prePulse(6).amps = 2.510;
    prePulse(6).analyzedPreResponse = 1;
    
    prePulse(7).elec = 20;
    prePulse(7).amps = 2.510;
    prePulse(7).analyzedPreResponse = 1;
    
    prePulse(8).elec = 21;
    prePulse(8).amps = 2.510;
    prePulse(8).analyzedPreResponse = 1;
    
    prePulse(9).elec = 22;
    prePulse(9).amps = [0.904 1.104 2.510];
    prePulse(9).analyzedPreResponse = [1 1 2];
    
    prePulse(10).elec = 23;
    prePulse(10).amps = 2.510;
    prePulse(10).analyzedPreResponse = 1;
    
    prePulse(11).elec = 24;
    prePulse(11).amps = 2.510;
    prePulse(11).analyzedPreResponse = 1;
    
    prePulse(12).elec = 26;
    prePulse(12).amps = [1.104 1.305 2.510];
    prePulse(12).analyzedPreResponse = [1 1 2];
    
    prePulse(13).elec = 27;
    prePulse(13).amps = 2.510;
    prePulse(13).analyzedPreResponse = 1;
    
    prePulse(14).elec = 28;
    prePulse(14).amps = [2.209 2.409 2.510];
    prePulse(14).analyzedPreResponse = [1 1 1];
    
    prePulse(15).elec = 29;
    prePulse(15).amps = 2.510;
    prePulse(15).analyzedPreResponse = 1;
    
    prePulse(16).elec = 32;
    prePulse(16).amps = 2.510;
    prePulse(16).analyzedPreResponse = 1;
    
    prePulse(17).elec = 33;
    prePulse(17).amps = [1.807 2.008 2.510];
    prePulse(17).analyzedPreResponse = [1 1 2];
    
    prePulse(18).elec = 36;
    prePulse(18).amps = 2.510;
    prePulse(18).analyzedPreResponse = 1;
    
    prePulse(19).elec = 41;
    prePulse(19).amps = 2.510;
    prePulse(19).analyzedPreResponse = [1 1 2];
end


%prePulse = extractSpatioTempProbeResults(pathToElecResp, stimElec, neuronID, prePulse, offsets, 'plotErfFits', true);
%keyboard
%save([pathToElecResp 'prePulse.mat'], 'prePulse')

load([pathToElecResp 'prePulse.mat'])

%% each pre elec plotted separately

nPreElecs = length(prePulse);

for ii = 1:nPreElecs
    plotSpatioTempProbeResults_1elec(prePulse, ii, offsets, pathToElecResp, neuronID, stimElec, ampRange, threshYLims)
end

%% summary figure - all preElecs plotted together

plotSpatioTempProbeResults_allElecs(prePulse, offsets, pathToElecResp, neuronID, stimElec, ampRange, threshYLims)