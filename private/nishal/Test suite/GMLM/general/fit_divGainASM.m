function fitASM = fit_divGainASM(X,Xtilde,Y,nSU,interval,f,fd)

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

initVal = 0.1*randn(idx,1);


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

finalVal = fminunc(@(y)fit_divgainASM_grad(y,X,Xtilde,Y,params,f,fd,interval),initVal,optim_struct);

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
