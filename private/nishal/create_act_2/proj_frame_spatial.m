function mov_new_frame=proj_frame_spatial(A,mov_fr,Ainv)

Filt_dim1=size(mov_fr,1);
Filt_dim2=size(mov_fr,2);

mov_fr=mov_fr(:);
%mov_null_old=mov_fr-A'*(A'\(A\(A*mov_fr))); % Doubt as its not full rank!!
mov_null=mov_fr-Ainv*A*mov_fr;
%norm(mov_null_old-mov_null)
%norm(A*mov_null)
mov_new_frame=reshape(mov_null,[Filt_dim1,Filt_dim2]);

end
