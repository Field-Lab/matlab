function res = rmse(d1,d2)

temp=single(d1(:)-d2(:)).^2;

res=sqrt( sum(temp(:))/prod(size(temp)) );