function [stas128x128,mov128x128]=preprocess128x128(stas,mov)
nSTAs=length(stas);
staLen = size(stas{1},4);
mov_len=size(mov,3)
dim1=size(mov,1);
dim2=size(mov,2);


% Mov
mov128x128=zeros(128,128,mov_len);
mov128x128(1:dim1,1:dim2,:)=mov;

stas128x128=cell(nSTAs,1);
for ista=1:nSTAs
stas128x128{ista}=zeros(128,128,1,staLen);
stas128x128{ista}(1:dim1,1:dim2,1,:)=stas{ista};
end


end