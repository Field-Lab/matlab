stimdim=4;
stimLen=100000;

% Gaussian
 sigma=1;
 x=randn(stimLen,stimdim)*sigma;

% symmetric Log cauchy
% meanx =0;
% sigma=1;
% scale_unif=2;
% scale_range=0.05;
% x=rand(stimLen,stimdim)*scale;
% x= exp(tan(pi*(x)/5));
% x=x-1;
% x=x.*sign(randn(size(x)));
% x=x*scale_range;
% figure;
% hist(x(:),100)

k=[1,0,0,1;
    0,1,1,0];

k_scale=1;

noSU = size(k,1);

lam=zeros(stimLen,1);
su_resp=zeros(stimLen,noSU);
for isu=1:noSU
    su_resp(:,isu) = exp(k_scale*k(isu,:)*x')';
    lam=lam+su_resp(:,isu);
end
N=@(x) max(x-3.2,0);
lam_n = N(lam);


S = poissrnd(lam_n);


soft_max_su =zeros(stimLen,noSU);
for isu=1:noSU
soft_max_su(:,isu) = su_resp(:,isu)./sum(su_resp,2);
end

sum(S)/(stimLen)



STx = x(S>0,:);
ST_soft_max = soft_max_su(S>0,:);
figure;
hist(ST_soft_max(:,1),20);

su1=soft_max_su(S>0,1)>0.5;
idx=1:4;
[NN,CC] = hist3(STx(su1,idx(logical(k(1,:)))),[20,20]);
figure;
contour(CC{1},CC{2},NN);
hold on;
plot([-5,5],[0,0],'r');
hold on;
plot([0,0],[-5,5],'b');

figure;
xxx=[-5:0.1:5];
[xx,hh]=hist(k_scale*k(1,:)*x',100);
plotyy(hh,xx/sum(xx),xxx,exp(xxx));

figure('Color','w');
scatter(su_resp(:,1),su_resp(:,2),20,1*ones(size(su_resp,1),1),'filled');
hold on
scatter(su_resp(S>0,1),su_resp(S>0,2),20,2*ones(size(su_resp(S>0,:),1),1),'filled');
xlabel('SU1');
ylabel('SU2');
title('Sub-unit inputs')


figure('Color','w');
scatter(x(:,1),x(:,4),20,1*ones(size(x,1),1),'filled');
hold on
scatter(x(S>0 & soft_max_su(:,1)>0.5,1),x(S>0 & soft_max_su(:,1)>0.5,4),20,2*ones(size(x(S>0 & soft_max_su(:,1)>0.5,:),1),1),'filled');
title('Pixel input')
xlabel('Pixel 1');
ylabel('Pixel 2');


%% GMLM 
[fitGLM,output] = fitGMLM_afterSTC_simplified(S,x',4,2) % works for sparse stimuli ..

%% spectral clustering of STx
npts = size(STx,1);
A=zeros(npts,npts);
sigma=1.5; % Adjust this!!

distance=zeros(npts,npts);

for i=1:npts
    for j=1:npts
    A(i,j)=exp(-norm(STx(i,:)-STx(j,:))^2 / (2*sigma^2));
    distance(i,j) = norm(STx(i,:)-STx(j,:));
    end
end

figure;
imagesc(A)

figure;
hist(A(:))

A(A<0.2)=0;
A=sparse(A);
figure;
spy(A)

d=sum(A,2);
D=diag(d);
L = eye(npts,npts) - (diag(sqrt(1./d)))*A*(diag(sqrt(1./d)));
% L =D-A;
L=sparse(L);

[E,lam]=eig(full(L));
lam=diag(lam);
[lam_sort,idx]=sort(lam,'ascend');

figure;
plot(lam_sort,'*');

E=E(:,idx);

kf=2;
U=E(:,1:kf);
for ipt=1:npts
U(ipt,:)=U(ipt,:)/norm(U(ipt,:));
end


[label,C]=kmeans(U,kf);
figure;
plot(U(label==1,1),U(label==1,2),'r*');
hold on
plot(U(label==2,1),U(label==2,2),'b+');

k_est=zeros(kf,stimdim);

for ik=1:kf
k_est(ik,:)= mean(STx(label==ik,:),1);
end

[label_sort,ii]=sort(label,'ascend');
A_new = A(ii,ii);
figure;
spy(A_new);

figure;
imagesc(A_new);

x2d = STx*k';
xd=x*k';
figure;
plot(xd(:,1),xd(:,2),'g.');
hold on;
plot(x2d(label==1,1),x2d(label==1,2),'r*');
hold on
plot(x2d(label==2,1),x2d(label==2,2),'b+');


figure;
icnt=0;
nc=stimdim;
X=x;
for dim1=1:stimdim
for dim2=1:stimdim
     icnt=icnt+1;
    if(dim1~=dim2)
    subplot(nc,nc,icnt);
    scatter(X(:,dim1),X(:,dim2),1.2,'g');
    hold on;
    scatter(STx(:,dim1),STx(:,dim2),1.2,'r');
    title(sprintf('dim1 %d dim2 %d',dim1,dim2));
    axis square
    end
end
end
