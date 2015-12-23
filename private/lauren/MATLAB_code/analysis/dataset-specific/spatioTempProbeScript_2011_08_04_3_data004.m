clear all

neuronID = 903;
stimElec = 1;
offsets = [0 10 20 40 80 160];
pathToElecResp = '/snle/lab/Experiments/Array/Analysis/2011-08-04-3/data004/';
ampRange = [0.2 0.8]; %used for plotting thresholds that are out of range and can't be determined
threshYLims = [0.2 0.9];


for x = 1; %for code folding purposes only
    prePulse(1).elec = 1;
    prePulse(1).amps = [0.301 0.402 2.008];
    prePulse(1).analyzedPreResponse = [1 1 2]; %1: analyzed, 2: not analyzed because almost all successes, 0: not analyzed for other reason
    
    prePulse(2).elec = 2;
    prePulse(2).amps = 2.008;
    prePulse(2).analyzedPreResponse = 1;
    
    prePulse(3).elec = 3;
    prePulse(3).amps = [1.807 2.008];
    prePulse(3).analyzedPreResponse = [1 1];
    
    prePulse(4).elec = 4;
    prePulse(4).amps = [0.703 0.803 2.008];
    prePulse(4).analyzedPreResponse = [1 1 2];
    
    prePulse(5).elec = 5;
    prePulse(5).amps = 2.008;
    prePulse(5).analyzedPreResponse = 1;
    
    prePulse(6).elec = 6;
    prePulse(6).amps = [1.606 1.807 2.008];
    prePulse(6).analyzedPreResponse = [1 1 1];
    
    prePulse(7).elec = 7;
    prePulse(7).amps = 2.008;
    prePulse(7).analyzedPreResponse = 1;
    
    prePulse(8).elec = 11;
    prePulse(8).amps = 2.008;
    prePulse(8).analyzedPreResponse = 1;
    
    prePulse(9).elec = 12;
    prePulse(9).amps = 2.008;
    prePulse(9).analyzedPreResponse = 1;
    
    prePulse(10).elec = 15;
    prePulse(10).amps = 2.008;
    prePulse(10).analyzedPreResponse = 1;
    
    prePulse(11).elec = 41;
    prePulse(11).amps = [1.506 1.707 2.008];
    prePulse(11).analyzedPreResponse = [1 1 1];
    
    prePulse(12).elec = 47;
    prePulse(12).amps = 2.008;
    prePulse(12).analyzedPreResponse = 1;
    
    prePulse(13).elec = 54;
    prePulse(13).amps = 2.008;
    prePulse(13).analyzedPreResponse = 1;
    
    prePulse(14).elec = 56;
    prePulse(14).amps = 2.008;
    prePulse(14).analyzedPreResponse = 1;
    
    prePulse(15).elec = 58;
    prePulse(15).amps = [0.703 2.008];
    prePulse(15).analyzedPreResponse = [1 2];
    
    prePulse(16).elec = 60;
    prePulse(16).amps = [1.004 1.205 2.008];
    prePulse(16).analyzedPreResponse = [1 1 2];
    
    prePulse(17).elec = 61;
    prePulse(17).amps = 2.008;
    prePulse(17).analyzedPreResponse = 1;
    
    prePulse(18).elec = 63;
    prePulse(18).amps = 2.008;
    prePulse(18).analyzedPreResponse = 1;
    
    prePulse(19).elec = 64;
    prePulse(19).amps = [0.904 1.104 2.008];
    prePulse(19).analyzedPreResponse = [1 1 2];
end


%prePulse = extractSpatioTempProbeResults(pathToElecResp, stimElec, neuronID, prePulse, offsets, 'plotErfFits', false);
%save([pathToElecResp 'prePulse.mat'], 'prePulse')

load([pathToElecResp 'prePulse.mat'])

%% each pre elec plotted separately

nPreElecs = length(prePulse);

for ii = 1:nPreElecs
    plotSpatioTempProbeResults_1elec(prePulse, ii, offsets, pathToElecResp, neuronID, stimElec, ampRange, threshYLims)
end

%% summary figure - all preElecs plotted together

plotSpatioTempProbeResults_allElecs(prePulse, offsets, pathToElecResp, neuronID, stimElec, ampRange, threshYLims)