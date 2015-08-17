function [LL,grad]=fit_divgainASM2_grad(y,X,Xtilde,Y,params,f,fd,interval)
nSU = length(params);

for isu = 1:nSU
params{isu}.k.value = y(params{isu}.k.idx);
params{isu}.a.value = y(params{isu}.a.idx);
params{isu}.b.value = y(params{isu}.b.idx);
params{isu}.c.value = y(params{isu}.c.idx);
params{isu}.sig.value = y(params{isu}.sig.idx);
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
    fdd = fd(params{isu}.k.value'*X+params{isu}.b.value)./(params{isu}.a.value'*Xtilde + params{isu}.sig.value);
    params{isu}.k.grad  = (X.*repmat(fdd,[size(X,1),1]))*differ;
    params{isu}.b.grad =fdd*differ;
    
    fda = -f(params{isu}.k.value'*X + params{isu}.b.value)./(params{isu}.a.value'*Xtilde + params{isu}.sig.value).^2;
    params{isu}.a.grad  = (Xtilde.*repmat(fda,[size(Xtilde,1),1]))*differ; 
    params{isu}.sig.grad = fda*differ;
    params{isu}.c.grad = sum(differ);
end

grad = gpuArray(0*y);
for isu=1:nSU
 grad(params{isu}.k.idx) = params{isu}.k.grad;
 grad(params{isu}.a.idx) = params{isu}.a.grad;
 grad(params{isu}.b.idx) = params{isu}.b.grad;
 grad(params{isu}.sig.idx) = params{isu}.sig.grad;
 grad(params{isu}.c.idx) = params{isu}.c.grad;
end

LL = gather(LL);
grad =gather(grad);
end