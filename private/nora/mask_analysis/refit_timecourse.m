function error = timecourse_error_function(timecourse; firing_rate, model_params, firing_rate)
fittedGLM_test = fittedGLM;
stimfilter           = fittedGLM.rawfit.spatialfilter * (timecourse);
stimfilter           = reshape(stimfilter, [ROI_length,ROI_length,length(paramind.X)]);
fittedGLM_test.linearfilters.Stimulus.Filter = stimfilter;
[gen_signal, ~] = glm_fit_NL_TEST(fittedGLM, fitspikes, fitmovie, 10);
predicted_firing = model_params(1)/(model_params(2)+exp(model_params(3)*gen_signal));
error = sum((predicted_firing - firing_rate).^2); % MSE
end

% fittedGLM.linearfilters.Stimulus.time_rk1'