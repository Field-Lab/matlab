base_filename = '/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk1_MU_PS_noCP_SUexp_init_p8IDp8fit/standardparams/WN_mapPRJ/2012-08-09-3/randomstart_';
SU_filters = zeros(8,9);
pooling_filter = zeros(8,169);
time_filter = zeros(8,30);
mu = zeros(8,1);
PS = zeros(8,120);
BPS = zeros(8,1);
for i = 1:8
   load([base_filename num2str(i) '/ONPar_841.mat'])
   SU_filters(i,:) = fittedGLM.SU_filter(:);
   pooling_filter(i,:) = fittedGLM.linearfilters.Stimulus.space_rk1(:);
   time_filter(i,:) = fittedGLM.linearfilters.Stimulus.time_rk1;
   mu(i,:) = fittedGLM.linearfilters.TonicDrive.Filter;
   PS(i,:) = fittedGLM.linearfilters.PostSpike.Filter;
   BPS(i,:) = fittedGLM.xvalperformance.logprob_glm_bpspike;
end