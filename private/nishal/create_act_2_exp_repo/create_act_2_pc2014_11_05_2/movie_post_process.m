function [mov_orig,mov_modify_new]=movie_post_process(mov_orig,mov_modify_new,mov_params)
if(strcmp(mov_params.mov_type,'bw')||strcmp(mov_params.mov_type,'bw-precomputed'))
scale=mov_params.deviation/quantile(abs(mov_modify_new(:)),1-mov_params.scaling_loss);
mov_modify_new=mov_modify_new*scale;


mov_modify_new(mov_modify_new>mov_params.deviation)=mov_params.deviation;
mov_modify_new(mov_modify_new<-mov_params.deviation) = -mov_params.deviation;


end

if(strcmp(mov_params.mov_type,'nsem'))
scale_plus=mov_params.deviation_plus/quantile(abs(mov_modify_new(:)),1-mov_params.scaling_loss);
scale_minus=-mov_params.deviation_minus/quantile(mov_modify_new(:),mov_params.scaling_loss);
scale=min([scale_plus,scale_minus]);

mov_modify_new=mov_modify_new*scale;


mov_modify_new(mov_modify_new>mov_params.deviation_plus)=mov_params.deviation_plus;
mov_modify_new(mov_modify_new<-mov_params.deviation_minus) = -mov_params.deviation_minus;
end


end