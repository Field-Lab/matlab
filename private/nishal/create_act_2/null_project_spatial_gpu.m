function [mov_orig,mov_new,stas_sp_current]=null_project_spatial_gpu(stas,mov,cell_params,matlab_cellids_datarun)
ncells=length(stas);
stas_sp_current=cell(length(stas),1);
% Choose best frame
if(cell_params.sta_spatial_method==1) % use 4th frame. Better if stas is clipped one. Best: cell_params.sta_spatial_method=1; usefits=2
    display('Finding best frame')
    for icell=1:length(stas)
        [V,I]=max(abs(squeeze(sum(sum(stas{icell}(:,:,1,:).^2,1),2))));
        I;
        %         plot((squeeze(sum(sum(stas{icell}(:,:,1,:),1),2))));
        stas_sp_current{icell}=stas{icell}(:,:,1,I);
    end
end
tic;
CellMasks_sp=cell_params.CellMasks;
if(cell_params.sta_spatial_method==4) %find temporal kernel by averaging and find ratio. Better if stas is clipped one. Best: cell_params.sta_spatial_method=4; usefits=2
    display('Find common waveform and regress');
    %     figure;
    parfor icell=1:length(stas)
        
        stas_sp_current{icell} = common_temporal_form(stas{icell},CellMasks_sp{icell})
        %
        % subplot(3,1,1);
        % imagesc(sum(stas{icell}(:,:,:,4),3));
        % axis image;
        % colormap gray
        % subplot(3,1,2);
        % imagesc(stas_sp_current{icell});
        % axis image;
        % colormap gray
        % title('spatial filter')
        % subplot(3,1,3);
        % plot(temporal_mean);
        % title('Temporal filter');
    end
end
toc;

% Choose rank 1 approximation. Best: sta_spatial_method=3 usefits = 0
tic;
if(cell_params.sta_spatial_method==3)
    display('Rank 1 approximation');
    for icell=1:ncells
        Amat=zeros(size(stas{icell},1)*size(stas{icell},2),size(stas{icell},4));
        for itime=1:size(stas{icell},4)
            x=stas{icell}(:,:,1,itime);
            Amat(:,itime)=x(:);
        end
        [u,s,v]=svds(Amat,1);
        Arecons=u(:,1)*s(1,1)*v(:,1)';
        stas_sp_current{icell}=reshape(Arecons(:,4),[size(stas{icell},1),size(stas{icell},2)]).*cell_params.CellMasks{icell};
        %
        % figure;
        % subplot(3,1,1);
        % imagesc(sum(stas{icell}(:,:,:,4),3));
        % axis image;
        % colormap gray
        % subplot(3,1,2);
        % imagesc(stas_sp_current{icell});
        % axis image;
        % colormap gray
        % title('spatial filter')
        % subplot(3,1,3);
        % plot(v(:,1));
        % title('Temporal filter');
    end
end
toc;


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

tic;
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

A = gpuArray(A);
[u,s,v]=svd(A,'econ');
Ainv=v*(s^-1)*u';

 mov=gpuArray(mov);
% Null each frame of movie
mov_new=0*mov;

for iframe=1:size(mov,3)
    
%         if(mod(iframe,100)==1)
%             iframe
%         end
    
    % mov_fr=mov(:,:,iframe);
    % mov_fr=mov_fr(:);
    % mov_null=mov_fr-A'*(A'\(A\(A*mov_fr)));
    % mov_new(:,:,iframe)=reshape(mov_null,[Filt_dim1,Filt_dim2]);
    mov_new(:,:,iframe)=proj_frame_spatial(A,mov(:,:,iframe),Ainv);
end


mov_orig=mov;
mov_orig = gather(mov_orig);mov_new = gather(mov_new);
%[proj,maxChange,averageChange]=computePixelChange(mov_new,mov_orig,A,Ainv,CellMasks_sp)
toc;


end