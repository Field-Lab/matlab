function [stas256x256,mov256x256]=preprocess256x256(stas,mov)
nSTAs=length(stas);
staLen = size(stas{1},4);
mov_len=size(mov,3)
dim1=size(mov,1);
dim2=size(mov,2);


% Mov
mov256x256=zeros(256,256,mov_len);
mov256x256(1:dim1,1:dim2,:)=mov;

stas256x256=cell(nSTAs,1);
for ista=1:nSTAs
stas256x256{ista}=zeros(256,256,1,staLen);
stas256x256{ista}(1:dim1,1:dim2,1,:)=stas{ista};
end


end