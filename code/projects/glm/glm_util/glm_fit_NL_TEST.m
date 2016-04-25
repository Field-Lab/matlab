%% NB 2015-05-01
function [gen_signal, firing_rate] = glm_fit_NL_TEST(fittedGLM, fitspikes, fitmovie, window)

gen_signal = glm_gen_signal_TEST(fittedGLM, fitmovie);
bins = length(gen_signal);
home_spbins  = ceil(fitspikes / fittedGLM.t_bin);
home_spbins = home_spbins(find(home_spbins < bins) );
firing = zeros(size(gen_signal));
firing(home_spbins) = 1;
firing_rate = conv(firing, gausswin(window), 'same')/sum(gausswin(window));

end
