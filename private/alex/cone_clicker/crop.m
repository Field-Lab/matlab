function crop(sbh)

global myPicture


myEl=round(ginput(2));

pdist(myEl)

figure
hist(double(myPicture(:)),50);

m=myPicture(myEl(2),myEl(1))



if 0
rect = round(getrect(sbh));

if rect(1)<1; rect(1)=1; end
if rect(2)<1; rect(2)=1; end
if (rect(1)+rect(3))>size(myPicture,2); rect(3)=size(myPicture,2)-rect(1); end
if (rect(2)+rect(4))>size(myPicture,1); rect(4)=size(myPicture,2)-rect(2); end

myPicture=myPicture(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3));

subplot(sbh)
colormap gray
imagesc(myPicture);

end