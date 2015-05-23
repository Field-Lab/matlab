cell='/OFFPar_346.mat';
exp='2013-10-10-0';

debug=[1 2 3 4 5 7 9 11 20];
NSEM_BPS=debug;
WN_BPS=debug;

type='NSEM_mapPRJ/';
[StimulusPars, ~, ~, ~] = Directories_Params_v23(exp, 'NSEM', 'PRJ');
debug_time_NSEM=debug*StimulusPars.slv.seconds_pernovelblock;

for i=1:length(debug)
    load(['/Volumes/Analysis/nora/NSEM/GLM_Output/debug_fixedSP_rk1_linear_MU_PS_noCP_p8IDp8/standardparams/' type exp '/debug_blocks_' num2str(debug(i)) cell]);
    NSEM_BPS(i)=fittedGLM.xvalperformance.glm_normedbits;
end

type='WN_mapPRJ/';
[StimulusPars, ~, ~, ~] = Directories_Params_v23(exp, 'WN', 'PRJ');
debug_time_WN=debug*StimulusPars.slv.seconds_pernovelblock;

for i=1:length(debug)
    load(['/Volumes/Analysis/nora/NSEM/GLM_Output/debug_fixedSP_rk1_linear_MU_PS_noCP_p8IDp8/standardparams/' type exp '/debug_blocks_' num2str(debug(i)) cell]);
    WN_BPS(i)=fittedGLM.xvalperformance.glm_normedbits;
end

hold on
plot(debug_time_NSEM,NSEM_BPS, 'b')
plot(debug_time_WN,WN_BPS, 'r')
legend ('NSEM', 'WN')
xlabel('Fitting Time (seconds)')
ylabel('BPS')

