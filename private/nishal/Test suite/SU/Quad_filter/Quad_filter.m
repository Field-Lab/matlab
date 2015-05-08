%% Generate data 
n=20000;
x1=normrnd(0,1,[n,1]);
x2=normrnd(0,1,[n,1]);

mn=zeros(4,1);
mn(1) = mean(x1);
mn(2) = mean(x2);

x1=x1-mn(1);
x2=x2-mn(2);

y1=x1.^2;
mn(3)=mean(y1);
y1=y1-mn(3);


y2=x2.^2;
mn(4) = mean(y2);
y2=y2-mn(4);

figure;
plot(x1,y1,'*');
x12 = x1.*x2;
mn(5)=mean(x12);
x12=x12-mn(5);
d=[x1,x2,y1,y2,x12];
C = d'*d;
[u,s,v]=svd(C);

M = (s^(-0.5))*u';
Minv = u*(s^(0.5));

d_whitened = M*d';
d_whitened = d_whitened';

d_whitened'*d_whitened

figure;
plot(d_whitened(:,2),d_whitened(:,1),'r*'); axis equal

%% 

dt=1/120;
%f=@(x)(1+2*x(:,1)).^2 + (1+2*x(:,2)).^2;
f=@(x)(2+2*x(:,1)+2*x(:,2)).^2;

%N=@(x) double((x-mean(x))>0).*(x-mean(x));
%N=@(x) 0.1*exp(0.3*x);
N=@(x) max(x-2,0);
%N=@(x) x.^2;

data_short = [x1,x2];

lambda = N(f(data_short))*dt;
resp = poissrnd(lambda);
resp(resp>0)=1;


sta= mean(data_short(resp==1,:),1)
sta_ex_non_white = mean(d(resp==1,:),1)
sta_whitened_data = mean(d_whitened(resp==1,:),1)

coeffs = Minv*sta_whitened_data'
coeffs2 = coeffs+mn

figure;

plot(data_short(:,1),data_short(:,2),'b*');
hold on;
plot(data_short(resp==1,1),data_short(resp==1,2),'r*');

figure;
idx1=2;
idx2=5;
plot(d_whitened(:,idx1),d_whitened(:,idx2),'b*');
hold on;
plot(d_whitened(resp==1,idx1),d_whitened(resp==1,idx2),'r*');

%% 

f_assumed = @(x) double((x-8.94)>0).*(x-8.94);
c=4;
b=[8;4];

cvx_begin
variable A(2,2) symmetric
variables b(2,1) c
subject to 

lam=data_short(1,:)*b + c + data_short(1,:)*A*data_short(1,:)';
for i=2:n
lam(i)=data_short(i,:)*b + c + data_short(i,:)*A*data_short(i,:)';
end

%lam_exp=0.1*exp(0.3*lam);
lam_exp = max(lam-2,0);
lam>=0

maximize (-sum(lam_exp*dt) + sum(lam(1,resp==1)-2))

cvx_end


%% 

cvx_begin
variable A(2,2) symmetric
variables b(2,1) c
subject to 

lam=data_short(1,:)*b + c + data_short(1,:)*A*data_short(1,:)';
for i=2:n
lam(i)=data_short(i,:)*b + c + data_short(i,:)*A*data_short(i,:)';
end

lam_exp=0.1*exp(0.3*lam);

%maximize (-sum(lam_exp*dt) + sum(lam(1,resp==1)))
maximize (sum(lam(1,resp==1)))

cvx_end



















