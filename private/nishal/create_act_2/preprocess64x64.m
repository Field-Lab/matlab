function [stas64x64,mov64x64]=preprocess64x64(stas,mov)
nSTAs=length(stas);
staLen = size(stas{1},4);
mov_len=size(mov,3)
dim1=size(mov,1);
dim2=size(mov,2);


% Mov
mov64x64=zeros(64,64,mov_len);
mov64x64(1:dim1,1:dim2,:)=mov;

stas64x64=cell(nSTAs,1);
for ista=1:nSTAs
stas64x64{ista}=zeros(64,64,1,staLen);
stas64x64{ista}(1:dim1,1:dim2,1,:)=stas{ista};
end


end