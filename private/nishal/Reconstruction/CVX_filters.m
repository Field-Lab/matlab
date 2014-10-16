addpath(genpath('/Volumes/Analysis/nishal/cvx'));

cvx_setup

%% 

filter_len=length(q);
lam=100;

cvx_begin 
variable filter_sp(filter_len)

minimize ((A*filter_sp-b)'*(A*filter_sp-b) + lam*norm(filter_sp,1));
cvx_end

%% 
% Faster Algorithm
%  normalize=0.0001*sum(A(:));
%  A=A/normalize;
%  b=b/normalize;

filter_len=size(A,2);
gam=0;
x=randn(filter_len,1);
a=1;
iter=0;
while a>10^-4
   iter=iter+1
   lk=0.5;%(1/iter);
   z=x-(lk*A'*(A*x-b));
   x_new=((z-lk*gam).*(z-lk*gam>0) - (-z-lk*gam).*(-z-lk*gam>0));
   a=norm(x-x_new)
   aa=norm(x)
   aaa=norm(A*x-b) 
    x=x_new;
end



%% L1/L2 norm ? 
filter_len=length(q);
cellFilter_len=(filter_len-1)/nCells;
lam=1000;
cvx_begin
variable filter_sp(filter_len)
subject to
f=0;
for iF=1:noCells
f=norm(filter_sp((iF-1)*cellFilter_len+2:iF*cellFilter_len+1),2);
end

minimize ((A*filter_sp-b)'*(A*filter_sp-b) + lam*f)
cvx_end

%%


filter_len=size(A,2);
gam=10;
x=10*randn(filter_len,1);
a=1;
iter=0;
while a>10^-4
   iter=iter+1
   lk=(1/iter);
   z=x-(lk*A'*(A*x-b));
   
   x_new=0*z;
   x_new(1)=z(1);
for iF=1:noCells
    z_cut=z((iF-1)*cellFilter_len+2:iF*cellFilter_len+1);
f=norm(z_cut,2);
x_new((iF-1)*cellFilter_len+2:iF*cellFilter_len+1)=(1-lk*gam/norm(z_cut))*z_cut.*(z_cut>lk*gam);
end
   %x_new=(x-gam*lk).*(x-gam*lk>0) - (-x-gam*lk).*(-x-gam*lk>0);
   
   
   a=norm(x-x_new)
   aa=norm(x)
    
    x=x_new;
end


