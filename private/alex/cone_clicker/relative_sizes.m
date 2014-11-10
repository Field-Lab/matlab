myPicture=imread('/Users/alexth/Dropbox/data/2014-04-15-5/Image4.jpg');
myPicture=double(myPicture(:,:,1));


myPicture2=imread('/Users/alexth/Dropbox/data/2014-04-15-5/Image17.jpg');
myPicture2=double(myPicture2(:,:,1));

comb=zeros(1200,1600,3);
comb(:,:,1)=myPicture/max(myPicture(:));
comb(:,:,3)=myPicture2/(max(myPicture2(:)));

figure
imagesc(comb)



figure
colormap gray
imagesc(myPicture/max(myPicture(:)))



figure
colormap gray
imagesc(myPicture2/max(myPicture2(:)))

myPicture=imread('/Volumes/Untitled/2014-08-19_test/2014-08-19 17.50.03.jpg');
myPicture=double(myPicture(:,:,1));


myPicture2=imread('/Volumes/Untitled/2014-08-19_test/Image7.jpg');
myPicture2=double(myPicture2(:,:,1));

figure
colormap gray
imagesc(myPicture)

figure
colormap gray
myPicture2=fliplr(myPicture2);
imagesc(myPicture2)



a=ginput(6);
a=diff(a(1,:)) % size of 1 checker

