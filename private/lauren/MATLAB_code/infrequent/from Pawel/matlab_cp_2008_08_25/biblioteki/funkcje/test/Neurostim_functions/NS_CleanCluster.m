function [Ind,a]=NS_CleanCluster(PCA_Coeff,Threshold);

SCoeff=size(PCA_Coeff);
SCoeff2=SCoeff;
Centre=zeros(1,SCoeff(2));

Centre=mean(PCA_Coeff); %Centre of mass for the cluster
STD=std(PCA_Coeff);
c=[Centre' PCA_Coeff']';
p=pdist(c);
distance=p(1:SCoeff(1)+1); %distance between mean and every element of PCA_Coeff;

Cond=zeros(SCoeff);
for i=1:SCoeff(2)
    a=PCA_Coeff(:,i)-Centre(i);
    b=find(abs(a)<Threshold*STD(i));
    Cond(b,i)=1;
end
Cond;
Ind=find(sum(Cond')==SCoeff(2));
%a1=find(SCoeff2<STD)';
%Ind=distance;