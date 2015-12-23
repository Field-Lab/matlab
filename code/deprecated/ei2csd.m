function csd = ei2csd(ei, neighbor_struct)

for i = 1:size(ei,2)
    csd(:,i) = lazyhex_laplacian(ei(:,i), neighbor_struct);
end