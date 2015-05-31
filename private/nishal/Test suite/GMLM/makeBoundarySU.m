function xxBig = makeBoundarySU(xx,nSU)
dim1=size(xx,1);
dim2=size(xx,2);

xxBig_fill = zeros(dim1*10,dim2*10);
for ix=1:size(xxBig_fill,1)
    for iy=1:size(xxBig_fill,2)
    xxBig_fill(ix,iy)= xx(floor((ix-1)/10)+1,floor((iy-1)/10)+1);
    end
end

f=[0,-1,0;
    -1,4,-1;
    0,-1,0];

xxBig = imfilter(xxBig_fill,f);
xxBig = xxBig~=0;
end