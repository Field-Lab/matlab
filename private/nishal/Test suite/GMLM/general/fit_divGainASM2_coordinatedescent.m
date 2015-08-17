function fitASM = fit_divGainASM2_coordinatedescent(X,Xtilde,Y,nSU,interval,f,fd,initVals)

% we get f(k'x+b)/(a'xt + sig) + c .. 

inpdim = size(X,1);
avginpdim = size(Xtilde,1);
T = size(Y,2);

params=cell(nSU,1);
idx=0;
for isu = 1:nSU
params{isu}.k.idx2 = idx+1 : idx+inpdim;
idx=idx + inpdim;
params{isu}.b.idx2 = idx+1:idx+1;
idx = idx+1;
params{isu}.c.idx2 = idx+1:idx+1;
idx = idx+1;
end
nrlen = idx;


idx=0;
for isu = 1:nSU
params{isu}.a.idx2 = idx+1:idx+avginpdim;
idx = idx + avginpdim;
params{isu}.sig.idx2 = idx+1:idx+1;
idx=idx+1;
end
drlen = idx;
%% initialize
% nr
initValnr = 0.01*rand(nrlen,1);
sta = X*Y'/sum(Y);
for isu=1:nSU
    if(isempty(initVals))    
        initValnr(params{isu}.k.idx2) =  sta + 0.1*randn(length(params{isu}.k.idx2),1)/norm(sta); % Careful when more SU
    else
        initValnr(params{isu}.k.idx2)=initVals{isu};
    end
end

% dr
% solve LP to make a'*xtilde + sig >0 
% A = [Xtilde',ones(size(Xtilde,2),1)];
% 
% optim_struct = optimset(...
%    'display','iter');  %'iter-detailed',... 
% 
% lpsol = linprog([zeros(size(Xtilde,1),1);0],-A,zeros(size(Xtilde,2),1),[],[],[],[],[],optim_struct);
initValdr = zeros(drlen,1);
for isu=1:nSU
initValdr(params{isu}.a.idx2) = 0;%lpsol(1:end-1);
initValdr(params{isu}.sig.idx2) = 1;
end

% initialize params
for isu=1:nSU
params{isu}.k.value = initValnr(params{isu}.k.idx2);
params{isu}.a.value = initValdr(params{isu}.a.idx2);
params{isu}.b.value = initValnr(params{isu}.b.idx2);
params{isu}.sig.value = initValdr(params{isu}.sig.idx2);
params{isu}.c.value = initValnr(params{isu}.c.idx2);
end

%%
optim_struct = optimset(...
   'derivativecheck','off',...
   'diagnostics','off',...  % 
   'display','iter-detailed',...  %'iter-detailed',... 
   'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
   'GradObj','on',...
   'largescale','on',...
   'Hessian','off',...
    'TolX',1e-8); 





X=gpuArray(X);
Xtilde = gpuArray(Xtilde);
Y = gpuArray(Y);

nr = initValnr;
dr = initValdr;
fval_prev=Inf;
togo=true;
while togo 
% fit gain
 
% numerator
[finalVal,fval] = fminunc(@(y)fit_divgainASM2_grad_gain_onlynr(y,X,Xtilde,Y,params,f,fd,interval),nr,optim_struct);

y=finalVal;
for isu = 1:nSU
params{isu}.k.value = y(params{isu}.k.idx2);
params{isu}.b.value = y(params{isu}.b.idx2);
params{isu}.c.value = y(params{isu}.c.idx2);
end
nr=finalVal;


% denominator
[finalVal,fval] = fminunc(@(y)fit_divgainASM2_grad_gain_onlydr(y,X,Xtilde,Y,params,f,fd,interval),dr,optim_struct);

y=finalVal;
for isu = 1:nSU
params{isu}.a.value = y(params{isu}.a.idx2);
params{isu}.sig.value = y(params{isu}.sig.idx2);
end
dr=finalVal;


if(abs(fval-fval_prev)>10^-6)
   togo=true; 
   fval_prev=fval;
else
togo=false;    
end

end

fitASM.params = params;

end
