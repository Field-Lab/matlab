cell='/ONPar_6858.mat';
exp='2012-09-27-3';

debug=[1 3 5 10 20 40 50];
NSEM_BPS=debug;
WN_BPS=debug;

type='NSEM_mapPRJ/';
[StimulusPars, ~, ~, ~] = Directories_Params_v23(exp, 'NSEM', 'PRJ');
debug_time_NSEM=debug*StimulusPars.slv.seconds_pernovelblock;
for i=1:length(debug)
    load(['/Volumes/Analysis/nora/NSEM/GLM_Output/debug_rk1_MU_PS_CP_p8IDp8/standardparams/' type exp '/debug_blocks_' num2str(debug(i)) cell]);
    NSEM_BPS(i)=fittedGLM.xvalperformance.glm_normedbits;
    NSEM.TempFilter{i}=fittedGLM.linearfilters.Stimulus.time_rk1;
    NSEM.PSFilter{i}=fittedGLM.linearfilters.PostSpike.Filter;
end

type='WN_mapPRJ/';
[StimulusPars, ~, ~, ~] = Directories_Params_v23(exp, 'WN', 'PRJ');
debug_time_WN=debug*StimulusPars.slv.seconds_pernovelblock;

for i=1:length(debug)
    load(['/Volumes/Analysis/nora/NSEM/GLM_Output/debug_rk1_MU_PS_CP_p8IDp8/standardparams/' type exp '/debug_blocks_' num2str(debug(i)) cell]);
    WN_BPS(i)=fittedGLM.xvalperformance.glm_normedbits;
    WN.TempFilter{i}=fittedGLM.linearfilters.Stimulus.time_rk1;
    WN.PSFilter{i}=fittedGLM.linearfilters.PostSpike.Filter;
end

hold on
plot(debug_time_NSEM,NSEM_BPS, 'b')
plot(debug_time_WN,WN_BPS, 'r')
legend ('NSEM', 'WN')
xlabel('Fitting Time (seconds)')
ylabel('BPS')

figure;
hold on
colors=(1:length(debug))/length(debug);
for i=1:length(debug)
    plot(WN.PSFilter{i}, 'Color',[1 1 1]*(1-colors(i)),'LineWidth',5*colors(i))
end
title('WN PS Filter')

figure;
hold on
colors=(1:length(debug))/length(debug);
for i=1:length(debug)
    plot(WN.TempFilter{i}, 'Color',[1 1 1]*(1-colors(i)),'LineWidth',5*colors(i))
end
title('WN Temporal Filter')

figure;
hold on
colors=(1:length(debug))/length(debug);
for i=1:length(debug)
    plot(NSEM.PSFilter{i}, 'Color',[1 1 1]*(1-colors(i)),'LineWidth',5*colors(i))
end
title('NSEM PS Filter')

figure;
hold on
colors=(1:length(debug))/length(debug);
for i=1:length(debug)
    plot(NSEM.TempFilter{i}, 'Color',[1 1 1]*(1-colors(i)),'LineWidth',5*colors(i))
end
title('NSEM Temporal Filter')