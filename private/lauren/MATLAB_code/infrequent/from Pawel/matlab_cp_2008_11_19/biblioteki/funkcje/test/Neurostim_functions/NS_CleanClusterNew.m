function [Ind,a]=NS_CleanCluster2(PCA_Coeff,Threshold);

SCoeff=size(PCA_Coeff)
Cond=zeros(SCoeff);
if SCoeff(1)>Threshold
    Cond=PCA_Coeff;
end