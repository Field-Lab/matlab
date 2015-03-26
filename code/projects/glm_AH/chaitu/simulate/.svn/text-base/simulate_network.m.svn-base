function [D hist_term input xout sta] = simulate_network(pars,boost,T,mu,sig)

K = zeros(pars.n*pars.Mk,pars.Nneurons);
ps = zeros(pars.Mhist,pars.Nneurons);
cp = zeros(pars.Mhist,pars.Nneurons*(pars.Nneurons-1));

for j=1:pars.Nneurons
    [kj ps(:,j) cp(:,(j-1)*(pars.Nneurons-1)+1:j*(pars.Nneurons-1))] = get_model_filters(pars,j);
    K(:,j) = reshape(kj',size(K,1),1);
end

% Simulate some data
[spikes hist_term input xout sta] = simulate_network_gwn_mex(pars.n,T,pars.dt,pars.b,mu,sig,boost.*K,ps,cp);
xout = reshape(xout,pars.n,T);

1;
% Construct a sparse binary vector with 1's at indices specified by spikes;
D = sparse(false(T,pars.Nneurons));
for j=1:pars.Nneurons
    D(spikes(spikes(:,2)==j,j)) = true;
end

fprintf('Recorded %d spikes.\n',size(spikes,1));