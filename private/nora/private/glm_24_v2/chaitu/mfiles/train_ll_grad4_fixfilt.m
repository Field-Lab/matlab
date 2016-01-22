% Function that computes the GRADIENT and HESSIAN (optional) of the log likelihood of the data D for the
% parameters p

% Arguments:

% basepars - basic parameters (see fitting_script.m)
% stimpars,trainpars - stimulus and training parameters, resp.
% cifs - the conditional intensity functions corresp. to the data wrt the
%        parameters p. This is the output of train_ll3.m

% Returns:
% ll_grad - gradient vector
% ll_hess - hessian matrix (optional)

function [ll_grad ll_hess] = train_ll_grad4_fixfilt(p,basepars,stimpars,trainpars,lcifs,cifs)

% Basic variables
%Ntrials = size(trainpars.D,2);
N = basepars.Nneurons;
Neff = size(trainpars.D,2);
ll_grad = zeros(size(p)); % this is the gradient vector

% Compute 'kernels' which are used as multipliers for computing
% gradients/hessian values. This is the factor that depends on the trial
% responses/cifs. It is multiplied by gradients/hessians of the input term
% (before the nonlinearity is applied).

%kernel = full((1.*trainpars.D) - (1.* ~trainpars.D).*((trainpars.dt.*cifs)./(1-(trainpars.dt.*cifs)))); % time x neuron
1;
fprime = basepars.Nprime(lcifs);
kernel = fprime./cifs.*full(double(trainpars.D(:,trainpars.baseneuron_idx))) - trainpars.dt.*fprime; % approx. for dt small
%kernel = full((double(trainpars.D(:,trainpars.baseneuron_idx))) - trainpars.dt.*cifs); % approx. for dt small
%kernel2 = full(-(1.* ~trainpars.D).*(trainpars.dt .* cifs)./((1-trainpars.dt.*cifs).^2)); % time x neuron
if (nargout > 1)
    %kernel2 = full(-trainpars.dt.*cifs); % approx for dt small
    fdprime = basepars.Ndoubleprime(lcifs);
    kernel2 = (fdprime.*cifs - fprime.^2)./(cifs.^2).*full(double(trainpars.D(:,trainpars.baseneuron_idx))) - fdprime.*trainpars.dt;
    ll_hess = zeros(length(p),length(p)); % this is the Hessian matrix
else
    kernel2 = [];
end
%kernel2_vec = squeeze(sum(kernel2,2)); % trial-summed version of the kernel: time x neuron


numpars = get_npars(basepars,Neff);        

T = size(trainpars.D,1); % time dim of the spike train (fine resolution!)
stimT = size(stimpars.x,2); % "" (coarse resolution)


for j=1:N
    trainpars.lgrad{j} = trainpars.kx(:,j)';
end

for n=1:N % Compute gradient/hessian of parameters corresponding to each neuron (all cross terms in Hessian are 0)

    idx = ((n-1)*numpars+1):(n*numpars); % indices for this neuron in the param vector p


    % Make a numpars x pars.maxt matrix where the jth column is the
    % gradient of the input term at time j for this trial (for the nth neuron)

    % Compute the gradients of the input terms wrt the filters: (2*(nspace+pars.Mk) x pars.maxt)
        
    % With multiple neurons, compute a matrix cpGrad where blocks of
    % basepars.nofilters_coupling rows are the gradients with respect to
    % the basis parameters for each neighboring neuron
    
    n_idx = trainpars.baseneuron_idx(n); % index of this neuron's spike train in trainpars.D
    neighbor_idx = [1:n_idx-1 n_idx+1:Neff]; % indices of the neighbor neuron's spike trains in trainpars.D
    1;
    if (Neff > 1) 
        cpGrad = zeros(basepars.nofilters_coupling*(Neff-1),T);
        counter = 1;
        for n2 = neighbor_idx % neuron n is the one being AFFECTED
%            n2rev = n2;
%            if (n2 > n)
%                n2rev = n2rev-1;
%            end
            cpGrad((counter-1)*basepars.nofilters_coupling+1:counter*basepars.nofilters_coupling,:) = trainpars.cpbasisGrad{n2}; % Add this neuron n2's spike train convolution
            counter = counter + 1;
        end
    end

    % Compute any extrinsic signal gradients, if necessary
    if (isfield(basepars,'ext_timepts') && ~isempty(basepars.ext_timepts))
        1;
       trainpars.extgrad{1} = sincinterpmtx(makeaxis(stimpars.dt,basepars.maxt),basepars.ext_timepts)';
    end    
    
    1;
    if (nargout == 1) % Evaluate the gradient only

        if (Neff > 1)
                [ll_grad(idx)] = multisample_gradHess_structs(basepars,n,kernel,[],trainpars,cpGrad);
        else
                [ll_grad(idx)] = multisample_gradHess_structs(basepars,n,kernel,[],trainpars);
        end
        continue;
    end
    
    1;
    if (nargout > 1) % Compute the the hessian if needed
        
        if (Neff > 1)
                [ll_grad(idx) ll_hess(idx,idx)] = multisample_gradHess_structs(basepars,n,kernel,kernel2,trainpars,cpGrad);
        else
                [ll_grad(idx) ll_hess(idx,idx)] = multisample_gradHess_structs(basepars,n,kernel,kernel2,trainpars);
        end
    end
end