 display('scaling the movie')
scale=mov_params.scale;
mov_modify_new=mov_modify_new*scale;
mov_orig = mov_orig*scale;


% Add means
mov_modify_new=mov_modify_new+mov_params.mean;
mov_orig=mov_orig+mov_params.mean;

% Clip
mov_modify_new(mov_modify_new>255)=255;
mov_modify_new(mov_modify_new<0) = 0;