function fitASM = fit_divGainASM2(X,Xtilde,Y,nSU,interval,f,fd)

% we get f(k'x+b)/(a'xt + sig) + c .. 

inpdim = size(X,1);
avginpdim = size(Xtilde,1);
T = size(Y,2);

params=cell(nSU,1);
idx=0;
for isu = 1:nSU
params{isu}.k.idx = idx+1 : idx+inpdim;

idx=idx + inpdim;

params{isu}.a.idx = idx+1:idx+avginpdim;
idx = idx + avginpdim;

params{isu}.b.idx = idx+1:idx+1;
idx = idx+1;

params{isu}.c.idx = idx+1:idx+1;
idx = idx+1;

params{isu}.sig.idx = idx+1:idx+1;
idx=idx+1;
end

% initialize
initVal = 0.1*rand(idx,1);
initVal(params{1}.k.idx) = X*Y'/sum(Y); % Careful when more SU
% solve LP to make a'*xtilde + sig >0 
A = [Xtilde',ones(size(Xtilde,2),1)];

optim_struct = optimset(...
   'display','iter');  %'iter-detailed',... 

lpsol = linprog([zeros(size(Xtilde,1),1);0],-A,zeros(size(Xtilde,2),1),[],[],[],[],[],optim_struct);
initVal(params{1}.a.idx) = lpsol(1:end-1);
initVal(params{1}.sig.idx) = lpsol(end);

optim_struct = optimset(...
   'derivativecheck','off',...
   'diagnostics','off',...  % 
   'display','iter-detailed',...  %'iter-detailed',... 
   'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
   'GradObj','on',...
   'largescale','on',...
   'Hessian','off',...
    'TolX',1e-8); 

%NOTE : Change Hessian to 'on' if we plan to use it .. 
X=gpuArray(X);
Xtilde = gpuArray(Xtilde);
Y = gpuArray(Y);

finalVal = fminunc(@(y)fit_divgainASM2_grad(y,X,Xtilde,Y,params,f,fd,interval),initVal,optim_struct);

y=finalVal;
for isu = 1:nSU
params{isu}.k.value = y(params{isu}.k.idx);
params{isu}.a.value = y(params{isu}.a.idx);
params{isu}.b.value = y(params{isu}.b.idx);
params{isu}.c.value = y(params{isu}.c.idx);
params{isu}.sig.value = y(params{isu}.sig.idx);
end

fitASM.params = params;

end
