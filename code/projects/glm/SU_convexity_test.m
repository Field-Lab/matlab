

p_SU = zeros(9,1);
test = 1.9:0.1:2.5;
val = zeros(length(test),1);

for i = 1:length(test)
   p_SU(5) = test(i);
   disp(i)
   val(i) = glm_SU_optimizationfunction_MSE(p_SU,SU_cov,pooling_weights,home_spbins,t_bin, non_stim_lcif);
end

plot(val)