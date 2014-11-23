function mov2 = mov3Dto4D(mov)

mov2=zeros(size(mov,1),size(mov,2),3,size(mov,3));

for itime=1:size(mov,3)
mov2(:,:,1,itime)=mov(:,:,itime);

mov2(:,:,2,itime)=mov(:,:,itime);

mov2(:,:,3,itime)=mov(:,:,itime);
end
end