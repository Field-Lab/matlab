m=reshape(k,400,400);
imagesc(m)
seq=[1 0 0 0 0 1 1 0 1 1 0 1 1 0 0 0 1 0 0 0 0 1 0 0];
seq=int2str(seq')'

find(k)
m=int2str(k')';
regexp(m,seq)

imagesc(sample_1_fine)

data=sample_1_fine(:,113:912,1);
data(data==60)=1;
colormap('gray')
imagesc(data)

seq=data(1:2:200);
seq=int2str(seq')'
regexp(m,seq)

imagesc(sample_1)
data=sample_1(:,113:912,1);
colormap('gray')
imagesc(data)
a=high_contr(1,:);
data(1:20:200)
seq=data(1:20:200);
seq=int2str(seq')'
m=int2str(a')';

imagesc(reshape(a,40,40))
figure
imagesc(data)




