function [mov_orig,mov_modify_new]=movie_post_process(mov_orig,mov_modify_new,mov_params)


if(strcmp(mov_params.mov_type,'bw')||strcmp(mov_params.mov_type,'bw-precomputed')||strcmp(mov_params.mov_type,'userProvided')||strcmp(mov_params.mov_type,'nsem'))
    
 if(~isfield(mov_params,'post_process_method'))   
 stretch_post_process_bw
 else

 if(strcmp(mov_params.post_process_method,'scale'))
 scale_post_process_bw
 else
 stretch_post_process_bw
 end
 
 end
 
end
% 
% if(strcmp(mov_params.mov_type,'nsem'))
% 
%  if(~isfield(mov_params,'post_process_method'))   
% stretch_post_process_nsem
%  else
% 
%  if(strcmp(mov_params.post_process_method,'scale'))
% scale_post_process_nsem % Not implemented yet!
%  else
% stretch_post_process_nsem
%  end
%  
%  end
% 
% 
% end


end