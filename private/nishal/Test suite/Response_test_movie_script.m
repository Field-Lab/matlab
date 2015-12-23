
figure;
% Linear response
cell_resp= Ax(stas,movieTest,movieL,1);

subplot(4,1,1);
plot(cell_resp);
title('Linear resp');

% Non linear response
NL_cell_resp=N(cell_resp);
subplot(4,1,2);
plot(NL_cell_resp);
title('Non Linear resp');

% Generate spikes @ movie frame rate resolution
spksTrial=zeros(movieL,nTrials);
for iTrial=1:nTrials
iTrial
spksTrial(:,iTrial)=poissrnd(NL_cell_resp);
end

subplot(4,1,3);
plot(spksTrial);
title('Binned Spikes');


% Or, Generate spikes @ other resolution
resolution=1/(120*10); % in seconds
lambdaCorrection = resolution/(1/120);
step=1/lambdaCorrection; % has to be integer
spksTrialsmall=zeros(movieL*step,nTrials);
for iTrial=1:nTrials
iTrial
icnt=0;
for itime=1:movieL
    for istep=1:step
        icnt=icnt+1;
        spksTrialsmall(icnt,iTrial)=poissrnd(NL_cell_resp(itime)*lambdaCorrection);
    end
end
end

spksTrialsmall(spksTrialsmall>0)=1; % remove all multiple spikes!

% plot rasters

figure;    
plotSpikeRaster(logical(spksTrialsmall'),'PlotType','vertline');
title('Raster');