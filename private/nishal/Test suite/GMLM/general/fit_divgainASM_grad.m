function [LL,grad]=fit_divgainASM_grad(y,X,Xtilde,Y,params,f,fd,interval)
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
  total_inp = total_inp + f((params{isu}.k.value'*X + params{isu}.b.value) ./(params{isu}.a.value'*Xtilde + params{isu}.sig.value) + params{isu}.c.value);    
end

LL = sum(total_inp)*interval/(120*T) - Y*log(total_inp')/T;

%% compute gradient
for isu = 1:nSU
    fder = (fd((params{isu}.k.value'*X + params{isu}.b.value) ./(params{isu}.a.value'*Xtilde + params{isu}.sig.value) + params{isu}.c.value).*(interval/(120*T)  - Y./(T*total_inp)))';
    params{isu}.k.grad  = (X./repmat(params{isu}.a.value'*Xtilde + params{isu}.sig.value,[size(X,1),1]))*fder;
    params{isu}.a.grad  = (Xtilde.*repmat((-(params{isu}.k.value'*X + params{isu}.b.value)./(params{isu}.a.value'*Xtilde + params{isu}.sig.value).^2),[size(Xtilde,1),1]))*fder; 
    params{isu}.b.grad = ((params{isu}.a.value'*Xtilde + params{isu}.sig.value).^(-1))*fder;
    params{isu}.sig.grad = (-(params{isu}.k.value'*X + params{isu}.b.value)./(params{isu}.a.value'*Xtilde + params{isu}.sig.value).^2)*fder;
    params{isu}.c.grad = sum(fder);
end

grad = 0*y;
for isu=1:nSU
 grad(params{isu}.k.idx) = params{isu}.k.grad;
 grad(params{isu}.a.idx) = params{isu}.a.grad;
 grad(params{isu}.b.idx) = params{isu}.b.grad;
 grad(params{isu}.sig.idx) = params{isu}.sig.grad;
 grad(params{isu}.c.idx) = params{isu}.c.grad;
end
end