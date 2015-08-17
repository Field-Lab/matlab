function [LL,grad]=fit_divgainASM2_grad_gain_onlydr(y,X,Xtilde,Y,params,f,fd,interval)
nSU = length(params);

for isu = 1:nSU
params{isu}.a.value = y(params{isu}.a.idx2);
params{isu}.sig.value = y(params{isu}.sig.idx2);
end


%% compute likelihood
T=length(Y);
total_inp =zeros(1,T);
for isu = 1:nSU
  total_inp = total_inp + f((params{isu}.k.value'*X + params{isu}.b.value)) ./(params{isu}.a.value'*Xtilde + params{isu}.sig.value) + params{isu}.c.value;    
end

total_inp = gather(total_inp);

LL = sum(total_inp)*interval/(120*T) - Y*log((total_inp'))/T;

%% compute gradient
for isu = 1:nSU
    differ =(interval/(120*T)  - Y./(T*total_inp))';

    fda = -f(params{isu}.k.value'*X + params{isu}.b.value)./(params{isu}.a.value'*Xtilde + params{isu}.sig.value).^2;
    params{isu}.a.grad  = (Xtilde.*repmat(fda,[size(Xtilde,1),1]))*differ; 
    params{isu}.sig.grad = fda*differ;
end

grad = gpuArray(0*y);
for isu=1:nSU
 grad(params{isu}.a.idx2) = params{isu}.a.grad;
 grad(params{isu}.sig.idx2) = params{isu}.sig.grad;
end

LL = gather(LL);
grad =gather(grad);
end