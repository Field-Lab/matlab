a=rand(1,1000000);
b=rand(1,1000000);
tic
z=fastcorr2(a,b,[10000 10000],16);
toc
a=rand(1,1000000);
b=rand(1,1000000);
tic
z1=fastcorr2c(a,b,[10000 10000],16);
toc