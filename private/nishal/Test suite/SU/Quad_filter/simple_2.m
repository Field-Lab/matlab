
stimLen=4;
b=[1,1,1,1]';
B = [1,0,1,0;
    0,1,0,0;
    0,0,0,1];
B=B'/sqrt(2); % use sqrt(2) , can recover 
A=B*B';

% If A too smaller than b, cant recover it .. some sort of transition on
% when I can recover A ?? ..

c=1;

n=20000;
data=randn(n,stimLen);

lam =data*b + c + sum((data*B).^2,2);
lam=lam';
lam_exp=10*exp(0.15*lam);


dt=1/120;
firingRate = lam_exp*dt;

Y = poissrnd(firingRate);
%Y(Y>0)=1;

sum(Y)/(n*dt)


%% CVX
data_short=data;
dt=1/120;
resp=Y;
n=length(resp);
gamma=0;
tic;
cvx_begin
variable An(3,3) semidefinite

variables bn(3,1) cn
subject to 

lam = data_short*bn + cn + sum((An*data_short').*data_short',1)';
lam=lam';
f=+(sum(10*exp(0.15*lam)*dt) - sum((0.15)*lam(resp==1))) + sum(gamma*abs(An(:)));

minimize (f)
cvx_end
toc;

%% 
% global rankB maskedMov binnedResponses masksz
% masksz = 3; 
% rankB=2;
% binnedResponses=Y';
% maskedMov=data';
%  nargs = 1 + masksz + masksz*rankB ;
%  x0(1)=c;
%  x0(2:2+masksz-1)=b;
%  x0(2+masksz:2+masksz+(masksz*rankB)-1) = B(:)';
% 
% %x0=ones(nargs,1)/10; % initialization very important !! 
% 
% 
%  options = optimoptions('fminunc','GradObj','on','Diagnostics','on','Display','iter-detailed');
% [x,fval] = fminunc(@quad_GLM_B,x0,options);
% 
% cn=x(1);
% bn=x(2:2+masksz-1);
% Bn=reshape(x(2+masksz:end),[masksz,rankB]);
% 
% 
%  data_out.A=Bn*Bn';
%  data_out.b=bn;
%  data_out.c=cn;

 %% 
 global maskedMov binnedResponses masksz gamma wt
 masksz = 4; 

binnedResponses=Y';
maskedMov=data';

 nargs = 1 + masksz + masksz^2 ;
 x0=randn(nargs,1)/10;
%x0=ones(nargs,1)/10;
%x0(1) = 0.4481;
%x0(2:2+masksz-1)=0.8009;
% 
%  x0(1)=c;
%  x0(2:2+masksz-1)=b;
%  A=B*B';
%  x0(2+masksz:2+masksz+(masksz*masksz)-1) = A(:)';
 
 
gamma=0;
wt=ones(masksz*masksz,1);
eps=0.1;
for iter=1:1
    
options = optimoptions('fminunc','GradObj','on','Diagnostics','on','Display','iter-detailed');
[x,fval] = fminunc(@quad_GLM_A,x0,options);
cn=x(1);
bn=x(2:2+masksz-1);
An=reshape(x(2+masksz:end),[masksz,masksz]);

wt = 1./(eps+x(2+masksz:end));
x0=x;
end


 data_out.A=An;
 data_out.b=bn;
 data_out.c=cn;
 
 %% make null stim
 n=20000; % lesser data in null ?
 
data=randn(n,stimLen);
weights = ones(1,stimLen);%data.b;
A_null=weights;
[u,s,v]=svd(A_null,'econ');
Ainv=v*(s^-1)*u';

mov_modify=0*data;
tic;
for iframe=1:n
    if(rem(iframe,1000)==1)
    iframe
    end
   mov_fr = data(iframe,:)';
mov_null=mov_fr-Ainv*A_null*mov_fr;
mov_modify(iframe,:)=mov_null';
end
toc;
data_orig=data;
data=mov_modify;
 
lam =data*b + c + sum((data*B).^2,2);
lam=lam';
lam_exp=10*exp(0.15*lam);


dt=1/120;
firingRate = lam_exp*dt;

Y = poissrnd(firingRate);
%Y(Y>0)=1;

sum(Y)/(n*dt)

%% 

 %global maskedMov binnedResponses masksz gamma wt 
 global  b_ratios
 b_ratios=ones(stimLen,1);
 masksz = 4; 
 binnedResponses=Y';
 maskedMov=data';

 nargs = 1 + 1 + masksz^2 ;
x0=randn(nargs,1)/10;
A0=reshape(x0(3:end),[masksz,masksz]);
A0=(A0+A0')/2;
x0(3:end)=A0(:);

%x0=ones(nargs,1)/10;
%x0(1) = 0.4481;
%x0(2:2+masksz-1)=0.8009;

%   x0(1)=c;
%   x0(2)=1;
%   A=B*B';
%   x0(3:2+1+(masksz*masksz)-1) = A(:)';
%  
 
gamma=0;
wt=ones(masksz*masksz,1);
eps=0.01;
for iter=1:1
    
options = optimoptions('fminunc','GradObj','on','Diagnostics','on','Display','iter-detailed');
[x,fval] = fminunc(@quad_GLM_A_b_fixed,x0,options);
cn=x(1);
bn=x(2)*b_ratios;
An=reshape(x(3:end),[masksz,masksz]);

wt = 1./(eps+x(3:end));
x0=x;
end


 data_out.A=An;
 data_out.b=bn;
 data_out.c=cn;
