
function u_spatial = reshape_vector2(u,masked_frame,indexed_frame)
dim1 = size(indexed_frame,1);
dim2 = size(indexed_frame,2);
u_spatial= zeros(dim1,dim2);

for ielem=1:length(u)
[r,c]=find(masked_frame(ielem)==indexed_frame);
u_spatial(r,c) = u(ielem);
end


end
