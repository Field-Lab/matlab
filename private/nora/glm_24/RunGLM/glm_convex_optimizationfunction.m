% AKHeitman 2014-05-13  this works!!  

function  [f grad Hess log_cif]= glm_convex_optimizationfunction(linear_params,covariates,spikebins,bin_duration)

% Initialize
p = linear_params;
COV = covariates;
dt = bin_duration;
spt = spikebins;
COV = COV;


%display(sprintf('tonic drive %d', p(1))  )
% Evaluate the building Blocks
lcif = p' * COV;
cif  = exp(lcif);
% Evaluate the function
f_eval = sum( lcif(spt) ) - dt * sum(cif);

% Evaluate the gradient
g_eval = sum(COV(:,spt),2)  - dt * ( COV * (cif') );

% Evaluate the Hessian
hessbase = zeros(size(COV));
for i_vec = 1:size(COV,1)
    hessbase(i_vec,:) = sqrt(cif) .* COV(i_vec,:) ;
end
H_eval = -dt * (hessbase * hessbase');


% Will run through fmin_unc.    Swithc the signs.
f       = -f_eval;
grad    = -g_eval;
Hess    = -H_eval;
log_cif = lcif;


end