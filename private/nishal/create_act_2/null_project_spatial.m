function [mov_orig,mov_new]=null_project_spatial(stas,mov,cell_params,matlab_cellids_datarun)

stas_sp_current=cell(length(stas),1);
if(cell_params.sta_spatial_method==1) % use 4th frame. Better if stas is clipped one.
    display('Finding best frame')
    for icell=1:length(stas)
        [V,I]=max(abs(squeeze(sum(sum(stas{icell}(:,:,1,:),1),2))));
        I
        plot((squeeze(sum(sum(stas{icell}(:,:,1,:),1),2))));
        stas_sp_current{icell}=stas{icell}(:,:,1,I);
    end
end

if(cell_params.sta_spatial_method==2) % fits spatial STA. better if stas is raw.
tic;
%fit_info=cell(length(stas),1);
precomputed.matlab_cellids=[];
precomputed.stas_sp=[];


try
    precomputed=load(cell_params.sta_spatial);
catch
display('Could not load sta_spatial');    
end

newcomputed=zeros(length(stas_sp_current),1);
parfor icell=1:length(stas)
    icell
    %fit_info{icell}=fit_sta_sequence(stas{icell}, 'fit_temporal',false,'fit_center',true,'fit_surround',true,'verbose',false);
    %full_fit{icnt} = sta_fit_function(fit_info{issta}.initial_params);
    %stas_sp{icell} = make_Gaussian_two_d(fit_sta_sequence(stas{icell}, 'fit_temporal',false,'fit_center',true,'fit_surround',true,'verbose',false));
    [stas_sp_current{icell},newcomputed(icell)]=fit_spatial_sta_for_nulling(stas{icell},precomputed,matlab_cellids_datarun(icell),cell_params.CellMasks{icell});
end

for icell=1:length(stas)
if(newcomputed(icell)==1)
    precomputed.matlab_cellids=[precomputed.matlab_cellids,matlab_cellids_datarun(icell)];
    precomputed.stas_sp{length(precomputed.stas_sp)+1}=stas_sp_current{icell};
end
end

try
    matlab_cellids=precomputed.matlab_cellids;
    stas_sp=precomputed.stas_sp;
    save(cell_params.sta_spatial,'matlab_cellids','stas_sp');%'precomputed.matlab_cellids','precomputed.stas_sp');
catch
display('Could not save');    
end

toc;
end

% Make A
Filt_dim1=size(stas_sp_current{1},1);
Filt_dim2=size(stas_sp_current{1},2);
ncells=length(stas);

A=zeros(ncells,Filt_dim1*Filt_dim2);
for icell=1:length(stas_sp_current)
A(icell,:)=stas_sp_current{icell}(:)';
end

if(rank(A)<min(size(A)))
display('Spatial STA matrix not well conditioned. Check selected cells.')
end

[u,s,v]=svd(A,'econ');
Ainv=v*(s^-1)*u';

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
mov_new(:,:,iframe)=proj_frame_spatial(A,mov(:,:,iframe),Ainv);
end

mov_orig=mov;

end