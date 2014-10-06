function [Types,Coeffs,Inc]=NS_ClusterSignatures(DATA,Dimensions,NumberOfClusters);

c1=princomp(DATA');
Coeffs=c1(:,1:Dimensions);
size(Coeffs);
Y = pdist(Coeffs,'euclidean'); 
Z = linkage(Y,'ward'); 
Inc=inconsistent(Z);
Types = cluster(Z,'maxclust',NumberOfClusters);