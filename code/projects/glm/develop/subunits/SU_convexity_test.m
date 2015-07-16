

p_SU = pstar_SU;
test = -1:0.1:-0.7;
val = zeros(length(test),1);

for i = 1:length(test)
   p_SU(5) = test(i);
   disp(i)
   val(i) = glm_SU_optimizationfunction_MSE(p_SU,SU_cov,pooling_weights,home_spbins,t_bin, non_stim_lcif);
end

plot(val)