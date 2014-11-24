function [mov_orig,mov_modify_new]=movie_post_process(mov_orig,mov_modify_new,mov_params)

mov_modify_new_masked = mov_modify_new(logical(repmat(mov_params.totalMaskAccept,[1,1,size(mov_modify_new,3)])));

if(strcmp(mov_params.mov_type,'bw')||strcmp(mov_params.mov_type,'bw-precomputed'))
    
scale=mov_params.deviation/quantile(abs(mov_modify_new_masked(:)),1-mov_params.scaling_loss);
mov_modify_new=mov_modify_new*scale;

% match modes of original and modified, change original
[h,x]=hist(abs(mov_modify_new_masked(mov_modify_new_masked>20)),100);
[maxfraq,maxvalue]=max(h);
mode_new = x(maxvalue);

[h,x]=hist(abs(mov_orig(:)),100);
[maxfraq,maxvalue]=max(h);
mode_orig = x(maxvalue);

orig_scale=mode_new/mode_orig;
mov_orig=orig_scale*mov_orig;

% Add means
mov_modify_new=mov_modify_new+mov_params.mean;
mov_orig=mov_orig+mov_params.mean;

% Clip
mov_modify_new(mov_modify_new>255)=255;
mov_modify_new(mov_modify_new<0) = 0;


end

if(strcmp(mov_params.mov_type,'nsem'))
scale_plus=mov_params.deviation_plus/quantile(abs(mov_modify_new_masked(:)),1-mov_params.scaling_loss);
scale_minus=-mov_params.deviation_minus/quantile(mov_modify_new(:),mov_params.scaling_loss);
scale=min([scale_plus,scale_minus]);

mov_modify_new=mov_modify_new*scale;




mov_modify_new=mov_modify_new+mov_params.mean;
mov_orig=mov_orig+mov_params.mean;

mov_modify_new(mov_modify_new>255)=255;
mov_modify_new(mov_modify_new<0) = 0;

end


end