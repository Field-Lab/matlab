

% RE FORMAT spkCondColl

spkCondCollformat(4).spksColl=[];
for icond=1:4
    spksColl=zeros(movie_time,nTrials);
    
    for itrial=1:nTrials
    for ispk=1:length(spkCondColl(icond).spksColl{itrial})
       frameNo= floor(double(spkCondColl(icond).spksColl{itrial}(ispk))*120/20000)+1;
       spksColl(frameNo,itrial)=1;
    end
    end
    spkCondCollformat(icond).spksColl=logical(spksColl)';
end


%%
% PSTH 



psthBinSize=10;
psthSmoothen=5;
for icond=1:4
[timeLogData{icond},psthData{icond}]=  psth_calc(( spkCondCollformat(icond).spksColl),psthBinSize,'nonoverlap');
[timeLogModel{icond},psthModel{icond}]=  psth_calc(( spkCondCollModel(icond).spksColl),psthBinSize,'nonoverlap');
psthModel{icond}=conv(psthModel{icond},(1/psthSmoothen)*ones(psthSmoothen,1),'same');
psthData{icond}=conv(psthData{icond},(1/psthSmoothen)*ones(psthSmoothen,1),'same');
end

figure;
for icond=1:4
    subplot(4,1,icond);
    plot(timeLogData{icond},psthData{icond},'b');
    hold on
    plot(timeLogModel{icond},psthModel{icond},'r');
    ylim([0,1.2*max(psthData{1})])
    title(sprintf('PSTH var: Data %f: LNP Model: %f',var(psthData{icond}),var(psthModel{icond})));
    legend('Data','LNP Model');
   
end
