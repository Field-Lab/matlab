function [mov_orig,mov_new]=null_project_spatial(stas,mov)

stas_sp=cell(length(stas),1);
for icell=1:length(stas)
    stas_sp{icell}=stas{icell}(:,:,1,4);
end

tic;
%fit_info=cell(length(stas),1);
parfor icell=1:length(stas)
    icell
    %fit_info{icell}=fit_sta_sequence(stas{icell}, 'fit_temporal',false,'fit_center',true,'fit_surround',true,'verbose',false);
    %full_fit{icnt} = sta_fit_function(fit_info{issta}.initial_params);
    stas_sp{icell} = make_Gaussian_two_d(fit_sta_sequence(stas{icell}, 'fit_temporal',false,'fit_center',true,'fit_surround',true,'verbose',false));
end
toc;
% Make A
Filt_dim1=size(stas_sp{1},1);
Filt_dim2=size(stas_sp{2},2);
ncells=length(stas);

A=zeros(ncells,Filt_dim1*Filt_dim2);
for icell=1:length(stas_sp)
A(icell,:)=stas_sp{icell}(:)';
end

% Null each frame of movie
mov_new=0*mov;
parfor iframe=1:size(mov,3)
%     if(mod(iframe,100)==1)
%         iframe
%     end
    
% mov_fr=mov(:,:,iframe);
% mov_fr=mov_fr(:);
% mov_null=mov_fr-A'*(A'\(A\(A*mov_fr)));
% mov_new(:,:,iframe)=reshape(mov_null,[Filt_dim1,Filt_dim2]);
mov_new(:,:,iframe)=proj_frame_spatial(A,mov(:,:,iframe));
end

mov_orig=mov;

end