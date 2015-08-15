% REFIT

% Convex Hull...just feed directly into minimization routine.. fast.. don't
% need to provide gradient or hessian
function obj_val = rescale_stim(scalars, lcif_stim0, home_spbins,t_bin)

offset    = scalars(1);
lcif_stim = scalars(2) * lcif_stim0;

total_lcif = offset + lcif_stim;
cif        = exp(total_lcif);
obj_val    = -(sum( total_lcif(home_spbins) ) - t_bin * sum(cif));
end

