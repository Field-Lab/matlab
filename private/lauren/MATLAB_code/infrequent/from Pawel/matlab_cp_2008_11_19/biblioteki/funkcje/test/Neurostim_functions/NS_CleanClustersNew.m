function TypesNew=NS_CleanClustersNew(PCA_Coeffs,Types,Threshold);

number=1;
while find(Types==number)
    number=number+1;
end
number=number-1; %number of clusters

TypesNew=zeros(size(Types));

for i=1:number
    Ind1=find(Types==i);
    if length(Ind1)<Threshold
        Types(Ind1)=0;
    end
end

TypesNew=Types;