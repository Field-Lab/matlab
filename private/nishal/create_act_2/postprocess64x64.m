function     [mov_orig,mov_modify_new] =postprocess64x64(mov_orig64x64,mov_modify_new64x64,mov)

mov_orig = mov_orig64x64(1:size(mov,1),1:size(mov,2),:);

mov_modify_new=mov_modify_new64x64(1:size(mov,1),1:size(mov,2),:);

end