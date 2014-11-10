
tl=double(imread('/Volumes/Data/Auxiliary/2014-04-10-2/spot/Image1.jpg'));

figure
colormap gray
imagesc(-tl+max(tl(:)))
tmp=tl;
tmp(tmp<60)=60;
figure
colormap gray
imagesc(tmp)

figure
plot(sort(tl(:)))


a=round(ginput(2));

myCone=tl(a(1,2):a(2,2),a(1,1):a(2,1));
figure
imagesc(myCone)
colormap gray

nnorm=myCone.*myCone;
nnorm=sum(nnorm(:));

cormat=zeros(1200-16,1600-16);
for i=1:1200-16
    for j=1:1600-16
        tmp=tl(i:i+16,j:j+16).*myCone;
        cormat(i,j)=sum(tmp(:));        
    end
end

cormat=cormat./nnorm;
figure
imagesc(cormat)
colormap gray
figure
imagesc(cormat(a(1,2):a(2,2),a(1,1):a(2,1)))
colormap gray

subs=tl(1:400,1:400);

figure
colormap gray
imagesc(subs)




b=subs;
LEN = 3;
THETA = 1;
PSF = fspecial('motion', LEN, THETA);
wnr1 = deconvwnr(subs, PSF, 0);
figure
imshow(wnr1);



b(b<50)=0;

acc=zeros([size(b),3]);
acc(:,:,1)=b;
acc(:,:,2)=subs;
acc=acc/max(acc(:));
figure
imagesc(acc)

acc=zeros(30000,2);
clear tmp
distance=10;
for i=1:3000
    [c,d]=max(b(:));
    [k,m]=ind2sub(size(tl),d);
    b(max(1,k-distance):min(k+distance,size(tl,1)),max(1,m-distance):min(m+distance,size(tl,2)))=0;
    
    D=pdist2([k,m],acc(1:i,:));
    tmp=min(D(D>0));
    
    if min(tmp)>30
        acc(i,:)=[k,m];
    end
%     imagesc(b)
%     drawnow
end





im6=double(imread('/Volumes/Data/Auxiliary/2014-04-10-2/spot/Image6.jpg'));
im7=double(imread('/Volumes/Data/Auxiliary/2014-04-10-2/spot/Image7.jpg'));
im8=double(imread('/Volumes/Data/Auxiliary/2014-04-10-2/spot/Image8.jpg'));
im9=double(imread('/Volumes/Data/Auxiliary/2014-04-10-2/spot/Image9.jpg'));

figure
colormap gray
tmp=im6/max(im6(:));
imagesc(tmp)
figure
hist(tmp(:),255)

comb=zeros([size(im6),3]);
comb(:,:,1)=im6/max(im6(:));
comb(:,:,2)=im7/max(im7(:));
comb(:,:,3)=im8/max(im8(:));

figure
imagesc(comb)


background = imopen(tmp,strel('disk',15));

% Display the Background Approximation as a Surface
figure, surf(double(background(1:8:end,1:8:end))),zlim([0 255]);
set(gca,'ydir','reverse');
figure
I2 = tmp - background;
colormap gray
imagesc(I2)

I3 = imadjust(I2);
figure
colormap gray
imagesc(I3);
level = graythresh(I3);
bw = im2bw(I3,level);
bw = bwareaopen(bw, 50);
figure
colormap gray
imagesc(bw);






%% "neural network"

n1=tl; % first layer is just the picture

tmp=tl;
tmp(tmp<55)=55;

figure
imagesc(tmp(20:end-20,20:end-20))
colormap gray

n2=conv2(tmp,myCone,'same');
figure
imagesc(n2(20:end-20,20:end-20))
colormap gray

comb=zeros(1200,1600,3);
comb(:,:,1)=tmp/max(tmp(:));
comb(:,:,2)=n2/max(n2(:))/2;
figure
imagesc(comb)

myCone % sample cone;
w12=
n2=



