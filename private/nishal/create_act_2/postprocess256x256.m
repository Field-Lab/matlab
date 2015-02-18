function     [mov_orig,mov_modify_new] =postprocess256x256(mov_orig256x256,mov_modify_new256x256,mov)

mov_orig = mov_orig256x256(1:size(mov,1),1:size(mov,2),:);

mov_modify_new=mov_modify_new256x256(1:size(mov,1),1:size(mov,2),:);

end