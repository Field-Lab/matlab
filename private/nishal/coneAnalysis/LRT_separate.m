cellID=2;

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

threshold=-0.4; % minimum violations
Connectivity = LRT_mat>threshold;
Connectivity(logical(eye(noCones)))=0