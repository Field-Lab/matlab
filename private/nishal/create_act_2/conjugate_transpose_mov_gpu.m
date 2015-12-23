function mov_ft_conj = conjugate_transpose_mov_gpu(mov_ft)

mov_ft=gpuArray(mov_ft);
mov_ft_conj=zeros(size(mov_ft),'gpuArray');
%
xLen=size(mov_ft,1);
yLen = size(mov_ft,2);
tLen = size(mov_ft,3);

for ix = 1:xLen
    for iy=1:yLen
        for it=1:tLen
            ixz=ix-1;
            iyz=iy-1;
            itz=it-1;
            mov_ft_conj(ixz+1,iyz+1,itz+1)=conj(mov_ft(mod(xLen-ixz,xLen)+1,mod(yLen-iyz,yLen)+1,mod(tLen-itz,tLen)+1));
        end
    end
end
end