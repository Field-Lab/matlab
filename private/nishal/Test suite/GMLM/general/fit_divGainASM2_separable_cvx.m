function fitASM = fit_divGainASM2_separable_cvx(X,Xtilde,Y,nSU,interval,f,fd,fdd,mask)

% we get f(k'x+b)/(a'xt + sig) + c .. 

inpdim = size(X,1);
avginpdim = size(Xtilde,1);
T = size(Y,2);

params=cell(nSU,1);

for isu = 1:nSU
idx=0;
params{isu}.k.idx2 = idx+1 : idx+inpdim;
idx=idx + inpdim;
end
nrlen = idx;
nrlenpersu = nrlen;


for isu = 1:nSU
    idx=0;
params{isu}.a.idx2 = idx+1:idx+avginpdim;
idx = idx + avginpdim;
params{isu}.sig.idx2 = idx+1:idx+1;
idx=idx+1;
end
drlen = idx;
drlenpersu = drlen;
%% initialize
mu=0.01;

% nr
initValnr = 0.01*rand(nrlen,1);
sta = X*Y'/sum(Y);
for isu=1:nSU
initValnr(params{isu}.k.idx2) =   0.1*randn(length(params{isu}.k.idx2),1)/norm(sta); % Careful when more SU
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
params{isu}.sig.value = initValdr(params{isu}.sig.idx2);
end

%%
optim_struct = optimset(...
   'derivativecheck','off',...
   'diagnostics','off',...  % 
   'display','iter-detailed',...  %'iter-detailed',... 
   'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
   'GradObj','on',...
   'largescale','on',...
   'Hessian','on',...
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

T=length(Y);
total_inp =zeros(1,T);
for isu = 1:nSU
  total_inp = total_inp + f((params{isu}.k.value'*X )) ./(params{isu}.a.value'*Xtilde + params{isu}.sig.value) ;    
end
total_inp = total_inp+mu;
fval_prev=Inf;


% numerator
togonr = true;
while togonr

 for isu=1:nSU
     nr=zeros(nrlenpersu,1);
     nr(params{isu}.k.idx2) = params{isu}.k.value;
    
     derikb = fd(params{isu}.k.value'*X)./(params{isu}.a.value'*Xtilde+params{isu}.sig.value)./log(total_inp);
     
   staf=X.*repmat(derikb,[size(X,1),1])*Y'/T;
   dnr = (params{isu}.a.value'*Xtilde + params{isu}.sig.value);
   [finalVal,fval] = fminunc(@(y)fit_divgainASM2_grad_gain_onlynr_cvx(y,X,Y,params,isu,f,fd,fdd,interval,staf,dnr),nr,optim_struct); 

    y=finalVal;
    params{isu}.k.value = y(params{isu}.k.idx2);    
    
    figure;
    u=zeros(80,40); u(mask) = gather(staf(1:end-1))%params{isu}.k.value(1:end-1);
    imagesc(u);
    axis image
    colormap gray
    
 end
 
T=length(Y);
total_inp =zeros(1,T);
for isu = 1:nSU
  total_inp = total_inp + (f(params{isu}.k.value'*X ) ./(params{isu}.a.value'*Xtilde + params{isu}.sig.value));    
end
total_inp = total_inp+mu;

fval  =  sum(total_inp)*interval/(120*T) - Y*log((total_inp'))/T

if abs(fval-fval_prev)>10^-6
togonr=true;
else
togonr=false;
end

fval_prev= fval;
end




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
