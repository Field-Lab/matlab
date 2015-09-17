function [labels,h] = cluster_spect(A,nSU);

dim = size(A,1);

A=A/max(A(:));
d= sum(A,2);
d = d.^(-1/2);
Dinv=diag(d);


L = eye(dim,dim) - Dinv*A*Dinv;

[E,Lam] = eig(L);
lam=diag(Lam);
[lam_sort,idx]=sort(lam,'ascend');
E=E(:,idx);


for ipt=1:dim
E(ipt,:) = E(ipt,:)/norm(E(ipt,:));
end

h=figure;
plot(lam_sort,'*');
%k=input('Insert number of clusters');
k=nSU;

E=E(:,1:k);
[labels,C]=kmeans(E,k);

end