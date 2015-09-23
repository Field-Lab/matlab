function [resp1,resp2] = framework_null_expt_2015_03_09_3(cellID,GIT_link,WN_orig,interval)
clear spatial_method
spatial_method=cell(4,1);
spatial_method{1}.struct_orig=[];
spatial_method{1}.struct_null=[];

icell_list=0;
% % for cellID= [2747];%]%[3152,3331,3692,4726,4921]%

    icell_list=icell_list+1
%fittedGLM=glm_fit_from_WNrun({3152,3331,3365,3620,3637,3692,3901,3902,3903,3904,3916,4129,4246,4291,4726,4789,4921,5059,5177,5326,5581,6006,6076,6391,6541,6725,6812,6826,6829,6856,7188,7532,7533,7651,7652,7726}, '2014-11-05-2/data009_nps', 'RGB-10-2-0.48-11111-32x32', 900, '/Volumes/Analysis/nora/nishal_glmfits/15min_rank2');
%save(sprintf('/Volumes/Analysis/nora/nishal_glmfits/15min_rank2/%d.mat',cellID),'fittedGLM');

load(sprintf('%s%d.mat',GIT_link,cellID));

%% Replace fitted linear filter with STA - better filter?
% 
% sta_filter=fittedGLM.cellinfo.WN_STA;
% cell_params.STAlen=30;
% sta_filt{1}=sta_filter;
% [new_stas,totalMaskAccept,CellMasks]=clipSTAs(sta_filt,cell_params);
% 
% sta_filter=squeeze(sum(sta_filter,3));
% xcoords = fittedGLM.linearfilters.Stimulus.x_coord;
% ycoords = fittedGLM.linearfilters.Stimulus.y_coord;
% totalMaskAccept=totalMaskAccept(ycoords,xcoords);
% sta_filter = sta_filter(ycoords,xcoords,:).*repmat(totalMaskAccept,[1,1,30]);
% sta_filter(:,:,1:14)=0;
% fittedGLM.linearfilters.Stimulus.Filter = sta_filter*max(abs(fittedGLM.linearfilters.Stimulus.Filter(:)))/max(abs(sta_filter(:)));
h1= figure;

%subplot(1,2,1);
imagesc(repelem(fittedGLM.linearfilters.Stimulus.Filter(:,:,6),20,20));
colormap gray
axis image
caxis([min(fittedGLM.linearfilters.Stimulus.Filter(:)),max(fittedGLM.linearfilters.Stimulus.Filter(:))]);
set(gca,'xTick',[]);set(gca,'yTick',[]);
%pause(1)
% 
% subplot(1,2,2);
% ssta_dummy = zeros(size(fittedGLM.linearfilters.Stimulus.Filter,1)^2,30);
% for itime =1:30
%     xx=fittedGLM.linearfilters.Stimulus.Filter(:,:,itime);
%     ssta_dummy(:,itime)=xx(:);
% end
% 
% s=svd(ssta_dummy);
% subplot(1,2,2);
% plot(s,'*');
% title(sprintf('Cell %d',cellID));

print(h1,'-depsc',sprintf('/Volumes/Lab/Users/bhaishahster/Spatial_null/Figures_EJ/cell_%d_linear_filter.eps',cellID));

%% Test cell
WNtime=120*24;
WNmovie =(double(rand(size(WN_orig,1),size(WN_orig,2),WNtime)>0.5)-0.5);
x=GLM_predict(fittedGLM, WNmovie+0.5, 50);
%plotraster(x.rasters.glm_sim,fittedGLM,'labels',true,'raster_length',24,'start_time',0)
figure;
plotSpikeRaster(logical(x.rasters.glm_sim));

    %% Generate response to WN
    
    WNtime=120*15*60/interval;
    WNmovie =double(rand(size(WN_orig,1),size(WN_orig,2),WNtime)>0.5)-0.5;
    WNmovie = repelem(WNmovie,1,1,interval);
    x=GLM_predict(fittedGLM, WNmovie+0.5, 1);
    
    %% Calculate STA 
    % WNmovie made 4 dimensional
    WNmov4D = zeros(size(WN_orig,1),size(WN_orig,2),3,WNtime*interval);
    for iframe=1:size(WNmovie,3);
    WNmov4D(:,:,1,iframe)=WNmovie(:,:,iframe);
    WNmov4D(:,:,2,iframe)=WNmovie(:,:,iframe);
    WNmov4D(:,:,3,iframe)=WNmovie(:,:,iframe);
    end
    mov_params.mov=WNmov4D;
    
    sta_params.Filtlen=30;
    sta_params.useTrial=1;
    
    cell_params.binsPerFrame=10;
    
    response.spksGen=x.rasters.glm_sim;
    aa=repmat([1:WNtime],[10,1]);
    response.mov_frame_number=aa(:);
    
    response = calculate_sta_ts(mov_params,response,sta_params,cell_params);
    WNSTA = response.analyse.STA;
         figure
         for itime=1:sta_params.Filtlen
             itime
         imagesc(squeeze((WNSTA(:,:,itime)))');colormap gray
         caxis([min(WNSTA(:)),max(WNSTA(:))]);
         colorbar
         pause(0.1);
         end

    %% Generate null movie from STA calculated above ? 
    use_fit_list=[2,0,0,2,0];
    sta_spatial_method_list=[1,2,3,4,0];
    ispatial_method=4%1:5%1:5
    use_fit_var=use_fit_list(ispatial_method) % 2, 0,0,2
    sta_spatial_method_var=sta_spatial_method_list(ispatial_method)%1,2 ,3,4
    if(ispatial_method~=5) 
      null_movie_compute_ts_spatial
    else
      null_movie_compute_ts
    end
      testmovie_filename='/Volumes/Lab/Users/bhaishahster/Spatial_null/Figures_EJ/18.rawMovie';
      testmovie=get_raw_movie(testmovie_filename,size(movie_full,3),1);
      testmovie=permute(testmovie,[2 3 1]);
       
    %% Generate rasters
x=GLM_predict(fittedGLM, double(testmovie)/255, 29);

figure;
plotSpikeRaster(x.rasters.glm_sim>0,'PlotType','vertline'); 

% make bin resolution / 10
resphr = x.rasters.glm_sim;
subSample=10;
resp=zeros(size(resphr,1),size(resphr,2)/subSample);
for ibin = 1:size(resphr,2)/subSample
resp(:,ibin) = sum(resphr(:,(ibin-1)*subSample + 1:(ibin-1)*subSample + subSample),2);
end

resp1 = resp(:,1:size(resp,2)/2);
resp2 = resp(:,size(resp,2)/2 + 1 : end);

end
