function     [mov_orig,mov_modify_new] =postprocess128x128(mov_orig128x128,mov_modify_new128x128,mov)

mov_orig = mov_orig128x128(1:size(mov,1),1:size(mov,2),:);

mov_modify_new=mov_modify_new128x128(1:size(mov,1),1:size(mov,2),:);

end