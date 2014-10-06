function [PCA_CoeffsNew,TypesNew]=NS_CleanClusters(PCA_Coeffs,Types,Threshold);

number=1;
while find(Types==number)
    number=number+1;
end
number=number-1; %number of clusters

PCA_CoeffsNew=PCA_Coeffs;
TypesNew=zeros(size(Types));
for i=1:number
    Ind1=find(Types==i);
    if length(Ind1)>1
        [Ind11,a1]=NS_CleanClusterNew(PCA_Coeffs(Ind1,:),Threshold);
        TypesNew(Ind1(Ind11))=i;    
    else
        TypesNew(Ind1)=i;
    end
end