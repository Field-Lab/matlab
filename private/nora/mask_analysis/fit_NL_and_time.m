function [model, t_init, increments] = fit_NL_and_time(iterations, fittedGLM, fitspikes, fitmovie, t_init, NL_init)
tic
for iter = 1:iterations
    disp(['iteration ' num2str(iter) ': fitting NL'])
    [gen_signal, firing_rate] = glm_fit_NL_TEST(fittedGLM, fitspikes, fitmovie, 1);
    model = fitnlm(gen_signal, firing_rate, 'y~b1./(b2+exp(b3*x1+b4))', NL_init);
    NL_init = model.Coefficients.Estimate;
    gen_signal_spatial = glm_gen_signal_spatial(fittedGLM, fitmovie);
    disp(['iteration ' num2str(iter) ': fitting timecourse'])
    new_time = fminunc(@(time_course)timecourse_error_function(time_course, firing_rate, model, gen_signal_spatial, 0*fittedGLM.linearfilters.TonicDrive.Filter), t_init);
    t_init = new_time;
    clear new_time;
    increments.NL{iter} = NL_init;
    increments.time{iter} = t_init;
    toc
end
end