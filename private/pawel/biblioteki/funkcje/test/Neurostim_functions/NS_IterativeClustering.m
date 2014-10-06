function y=NS_IterativeClustering(Data,Dimensions,NumberOfClusters);

threshold=1;

NumberOfClusters=2;
[Types,PCA_Coeffs,Inc]=NS_ClusterSignatures(Data,Dimensions,NumberOfClusters);

Ind1=find(Types==1);
[Ind11,a1]=NS_CleanCluster(PCA_Coeffs(Ind1,:),threshold);
Ind2=find(Types==2);
[Ind12,a2]=NS_CleanCluster(PCA_Coeffs(Ind2,:),threshold);

STD=std(PCA_Coeff(Ind11,:));
STD=std(PCA_Coeff(Ind12,:));