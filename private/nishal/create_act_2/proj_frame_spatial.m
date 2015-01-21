function mov_new_frame=proj_frame_spatial(A,mov_fr)

Filt_dim1=size(mov_fr,1);
Filt_dim2=size(mov_fr,2);

mov_fr=mov_fr(:);
mov_null=mov_fr-A'*(A'\(A\(A*mov_fr)));
mov_new_frame=reshape(mov_null,[Filt_dim1,Filt_dim2]);

end
