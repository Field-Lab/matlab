function  [f, grad, Hess, log_cifv] = glm_SU_optim(linear_params, convex_cov, stimulus, spikebins, bin_duration)

%% Initialize
dt = bin_duration;
spt = spikebins;
n_bins = size(convex_cov);
n_params = length(p);
n_time = length(stimulus);
n_lin_params = size(convex_cov,1);

% Initialize
p = linear_params(1:n_lin_params);
COV = convex_cov;
dt = bin_duration;
spt = spikebins;

% Find Conditional Intensity and its log
lcif_nonstim = p' * COV;
idx = n_lin_params;
local_stim = conv2(stimulus, reshape(linear_params((end-8):end), [3,3]), 'same');
SU_drive = log(1-exp(reshape(local_stim, [size(stimulus,1)*size(stimulus,2), n_time])));
lcif_stim = linear_params(idx:(end-9))*SU_drive;
cif = exp(lcif_stim + lcif_nonstim);

% Evaluate the objective function (monotonic in log-likelihood)
f_eval = sum( lcif(spt) ) - dt * sum(cif);

% Evaluate the gradient
g_eval(1:n_lin_params) = sum(COV(:,spt),2)  - dt * ( COV * (cif') );
dlcif = zeros(n_params - n_lin_params, n_time);
for i= (n_lin_params+1):(n_params-9)
   dlcif(i-n_lin_params, :) =  SU_drive(i-n_lin_params,:);
   temp = imresize(dlcif(i-n_lin_params, :), 10);
   g_eval(i) = sum(dlcif(spt))- dt*sum(cif * dlcif);
end
for i=(n_params-8):n_params
   dlcif =  SU_drive(i-n_lin_params,:);
   g_eval(i) = sum(dlcif(spt))- dt*sum(cif * dlcif);
end

% Evaluate the hessian
hessbase_nonstim = zeros(size(COV));
for i_vec = 1:size(COV,1)
    hessbase_nonstim(i_vec,:) = sqrt(cif) .* COV(i_vec,:) ;
end

SU_idx = 1;
pool_idx = 1;
location = zeros(121, 9, 2);
for i_SU = 1:3
    for j_SU = 1:3 
        for i_pool = 1:11
            for j_pool = 1:11
                location(pool_idx, SU_idx,:) = [i_pool, j_pool] + [i_SU,j_SU] - [3,3];
                pool_idx = pool_idx + 1;
            end
        end
        SU_idx = SU_idx+1;
        pool_idx = 1;
    end
end

for i_SU = 1:9
    for j= (n_lin_params+1):(n_params-9)
        dlcif = stimulus(location(j-n_lin_params,i_SU, :),:)./(exp(-local_stim)+1);
    end
    for j_SU = 1:9
        dlcif = linear_params(idx:(end-9))*stimulus(location(:,i_SU),:).*stimulus(location(:,j_SU),:).*(exp(-local_stim)+1)^(-2).*(exp(-local_stim)+1);
    end
end

H_eval = -dt * (hessbase * hessbase');


% Switch signs because using a minimizer  fmin
f       = -f_eval+f_penalty;
grad    = -g_eval+del_penalty';
Hess    = -H_eval;
log_cif = lcif;
end
