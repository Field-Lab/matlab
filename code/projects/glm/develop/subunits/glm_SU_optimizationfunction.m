function [f, grad, Hess, lcif] = glm_SU_optimizationfunction(p,filtertype,paramind,glm_covariate_vec,X_frame,frame_shifts, bpf, home_spbins,t_bin)

% non stimulus stuff
[f_nonstim, grad_nonstim, Hess_nonstim, lcif_nonstim] = glm_convex_optimizationfunction(p, glm_covariate_vec, spikebins, bin_direction);

% stimulus stuff

end