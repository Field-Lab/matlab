function p = blobdetector(im, side, min_dist,option,th,map);

% function p = blobdetector(im, side, min_dist,option,th,map);
%
% Input
%-------------------------------------------------------------
% im = image to count cells
% side: filter size - choose the diameter (or slightly larger) of the blob (in pixels)
% min_dist : minimum distance between peaks
% option : 1 if the peak to detect is dark
%          0 if the peak to detect is bright (e.g. Topro stained image)
% map : mask of layer
% th : threshold for filter output
      %- suggested to choose -0.05 or 0

% Output
%-------------------------------------------------------------
% p: returns the number of peaks detected.
% resultImaage: under current directory as 
% Result(filter size,min_dist)_p_totalarea

% Works only with gray-scale images. If RGB image is input, the image will
% be converted to gray image

%
tic

% specified result path
resultPath = [pwd,'\'];


if nargin < 4
    option = 0;
end

if nargin < 5
    th = 0;
end

% if map is not specified, assume a whole image as a map
if nargin < 6
    map = ones(size(im,1),size(im,2));
end
if size(im,3) ==1
    img = im;
else
    img = rgb2gray(im);
end

if option==0
    img = 255-img;
end

% apply the filter
ac = lapofgau(img,side);
ac(find(ac<th))=th;
ac = ac-th;

scaling = 1;
ac2=ac;
figure; imagesc(ac); axis image; axis off

% find local maxima from filter output
[i,j,val] = find_local_max_2D(ac2,[],floor(side/4/scaling),inf,min_dist,[],[1 1],map);

area = length(find(map==255 |map==1 ));

p = 0;
ind=[];
tmp = [i' j' ];
if size(im,3) == 1 % gray
    result(:,:,1) = im;
    result(:,:,2) = im;
    result(:,:,3) = im;
else
    result = im;
end

for k = 1:length(i)
    p = p+1;
    ind = [ind ; tmp(k,:)];
    result(tmp(k,1), tmp(k,2),1)= 255;
    result(tmp(k,1), tmp(k,2),2)= 255;
    result(tmp(k,1), tmp(k,2),3)= 0;
end


toc
figure;
imagesc(im);axis image; colormap gray; axis off
hold on; plot(scaling*ind(:,2), scaling*ind(:,1), '.r');
title([num2str(p),' blobs (area:',num2str(area),')']);
imwrite(result,[resultPath,'Result(',num2str(side),',',num2str(min_dist),')_',num2str(p),'_',num2str(area),'.tif'],'tif');

