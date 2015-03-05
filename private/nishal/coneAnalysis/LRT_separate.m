cellID=9;
figure;
imagesc(datarun.stas.rfs{cellID});


cones= find(datarun.cones.weights(:,cellID)==1);
coneInp = datarun.cone_inputs(:,cones);

y = double(datarun.spike_rate(cellID,:));

figure;
plot(coneInp)

f=  @(x) 1./(1+exp(-5*x)); % @(x) exp(1.2*x)%1./(1+exp(-5*x));

figure;
[X,N]=hist(coneInp(:));
bar(N,X/sum(X));
hold on
plot([-0.5:0.1:0.5],f([-0.5:0.1:0.5]),'r');

noCones=length(cones);
LRT_mat=[];
for cone1=1:noCones
for cone2=1:noCones;
    
if(cone1==cone2)
continue;
end
x1 = coneInp(:,cone1);
x2 = coneInp(:,cone2);

pnull_1 = f((x1+x2)/2);
palt_1 = (f(x1)+f(x2))/2;

evi_null = y*log(pnull_1) + (1-y)*log(1-pnull_1);
evi_alt =  y*log(palt_1) + (1-y)*log(1-palt_1);

LRT = evi_null-evi_alt;

LRT_mat(cone1,cone2)=LRT;

end
end

threshold=-166; % choose threshold according to minimum violations
Connectivity = LRT_mat>threshold;
Connectivity(logical(eye(noCones)))=0
% Cones 2, 5 for cell ID 2 must be connected
%%

coneInp_SpkTrig = coneInp(logical(y'),:);
covar = double(coneInp_SpkTrig'*coneInp_SpkTrig)/sum(y>0);

idx=logical(1-eye(noCones));
lam=0.15;
cvx_begin
variable c_sp(noCones,noCones) semidefinite
subject to

minimize (norm(covar-c_sp)+lam*(sum(abs(c_sp(idx)))) ) 
cvx_end


Connectivity2 = (abs(c_sp)>0.0001);
Connectivity2(logical(eye(noCones)))=0


figure;
subplot(1,4,1);
covar(logical(eye(noCones)))=0;
imagesc(covar);
axis image
colorbar


subplot(1,4,2);
c_sp(logical(eye(noCones)))=0;
imagesc(c_sp);
axis image
colorbar


subplot(1,4,3);
Connectivity2(logical(eye(noCones)))=0;
imagesc(Connectivity2);
axis image
colorbar

subplot(1,4,4);
Connectivity(logical(eye(noCones)))=0;
imagesc(Connectivity);
axis image
colorbar
