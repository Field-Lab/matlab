% stretch_post_process_nsem

display('Stretching the movie')

mov_modify_new_masked = mov_modify_new(logical(repmat(mov_params.totalMaskAccept,[1,1,size(mov_modify_new,3)])));

scale_plus=mov_params.deviation_plus/quantile(abs(mov_modify_new_masked(:)),1-mov_params.scaling_loss);
scale_minus=-mov_params.deviation_minus/quantile(mov_modify_new(:),mov_params.scaling_loss);
scale=min([scale_plus,scale_minus]);

mov_modify_new=mov_modify_new*scale;




% match modes of original and modified, change original
[h1,x1]=hist(abs(mov_modify_new(abs(mov_modify_new)>5)),100);
[maxfraq,maxvalue]=max(h1);
mode_new = x1(maxvalue);

[h2,x2]=hist(abs(mov_orig(abs(mov_orig)>5)),100);
[maxfraq,maxvalue]=max(h2);
mode_orig = x2(maxvalue);

figure;
plotyy(x2,h2,x1,h1);

orig_scale=mode_new/mode_orig;
mov_orig=orig_scale*mov_orig;

mov_modify_new=mov_modify_new+mov_params.mean;
mov_orig=mov_orig+mov_params.mean;

mov_modify_new(mov_modify_new>255)=255;
mov_modify_new(mov_modify_new<0) = 0;
