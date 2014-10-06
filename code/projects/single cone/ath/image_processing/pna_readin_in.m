array_cones=double(imread('/Volumes/Data/Auxiliary/2014-04-15-5/alignment_after_and_fluo_red_PNA/Image4.jpg'));
figure
colormap gray
imagesc(array_cones)

edited=double(imread('/Volumes/Data/Auxiliary/2014-04-15-5/alignment_after_and_fluo_red_PNA/edited_Image4.jpg'));
edited=edited(:,:,1);
figure
colormap gray
imagesc(edited)

myArray=imresize(edited,1.453);
myArray=myArray(1:1200,1:1600);
myArray=myArray-min(myArray(:));
figure
colormap gray
imagesc(myArray)


pna_cones=double(imread('/Volumes/Data/Auxiliary/2014-04-15-5/alignment_after_and_fluo_red_PNA/Image12.jpg'));
figure
colormap gray
imagesc(pna_cones)


figure
merged_cones=zeros(1200,1600,3);
merged_cones(:,:,1)=array_cones/max(array_cones(:));
merged_cones(:,:,2)=pna_cones/max(pna_cones(:))/2;
merged_cones(:,:,3)=myArray/max(myArray(:));

imagesc(merged_cones)
