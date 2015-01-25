%% add path to nora's folder for GLM code
location_of_git_repo='/home/vision/Nishal/matlab/';
addpath(genpath([location_of_git_repo '/private/nora']));

% Java library
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');
addpath(genpath('../../create_act2'));
addpath(genpath('../GLM'));

%%

% clear spatial_method
% spatial_method=cell(4,1);
% spatial_method{1}.struct_orig=[];
% spatial_method{1}.struct_null=[];

%icell_list=0;
for cellID=[4246,4291,4726,4789,4921,5059,5177,5326,5581,6006,6076,6391,6541,6725,6812,6826,6829,6856,7188,7532,7533,7651,7652,7726];%[3152,3331,3365,3620,3637,3692,3901,3902,3903,3904,3916,4129]%
try
    icell_list=icell_list+1
%fittedGLM=glm_fit_from_WNrun({cellID}, '2014-11-05-2/data009_nps', 'RGB-10-2-0.48-11111-32x32', 900, '~/Nishal/GLM_fits');
load(sprintf('/Volumes/Analysis/nora/nishal_glmfits/15min/%d.mat',cellID));
%% Test cell
WNtime=120*24;
WNmovie =double(rand(32,32,WNtime)>0.5)-0.5;
x=GLM_predict(fittedGLM, WNmovie, 50);
plotraster(x,fittedGLM,'labels',true,'raster_length',24,'start_time',0)
figure;
plotSpikeRaster(logical(x.rasters.glm_sim))

    %% Generate response to WN
    WNtime=120*60*30;
    WNmovie =double(rand(32,32,WNtime)>0.5)-0.5;
    x=GLM_predict(fittedGLM, WNmovie, 1);
    
    %% Calculate STA 
    % WNmovie made 4 dimensional
    WNmov4D = zeros(32,32,3,WNtime);
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
         imagesc(squeeze((WNSTA(:,:,itime)))');colormap gray
         caxis([min(WNSTA(:)),max(WNSTA(:))]);
         colorbar
         pause(1/120)
         end

    %% Generate null movie from STA calculated above ? 
    use_fit_list=[2,0,0,2,0];
    sta_spatial_method_list=[1,2,3,4,0];
    for ispatial_method=1:5%1:5
    use_fit_var=use_fit_list(ispatial_method) % 2, 0,0,2
    sta_spatial_method_var=sta_spatial_method_list(ispatial_method)%1,2 ,3,4
    if(ispatial_method~=5) 
      null_movie_compute_ts_spatial
    else
      null_movie_compute_ts
    end
      testmovie_filename='~/Nishal/TS_data/18.rawMovie';
      testmovie=get_rawmovie(testmovie_filename,2880);
      testmovie=permute(testmovie,[2 3 1]);
       
    %% Generate rasters
x=GLM_predict(fittedGLM, testmovie, 100);

plotraster(x,fittedGLM,'labels',true,'raster_length',24,'start_time',0)
% figure;
% plotSpikeRaster(logical(x.rasters.glm_sim))
original_struct= sqrt(var(x.rate(1:floor(end/2))));
null_struct=sqrt(var(x.rate(floor(end/2)+1:end)));
title(sprintf('Original struct %f , Null struct %f',original_struct,null_struct));

spatial_method{ispatial_method}.struct_orig(icell_list)=original_struct;
spatial_method{ispatial_method}.struct_null(icell_list)=null_struct;
    end
    
catch
end


end


%% 
cols='rgbkm';
figure;
subplot(1,2,1)
for ispatial_method=1:5
loglog(spatial_method{ispatial_method}.struct_orig,spatial_method{ispatial_method}.struct_null,'*','Color',cols(ispatial_method));
hold on
end
loglog([0.1,max(spatial_method{ispatial_method}.struct_orig)],[0.1,max(spatial_method{ispatial_method}.struct_orig)]);
legend('best frame','Diff of Gaussian','Low rank approx','common temporal','spatio-temporal','45 degree line','Location','best');
xlabel('Structure in Original');
ylabel('Sructure in Null');

subplot(1,2,2)

for ispatial_method=1:5
plot(spatial_method{ispatial_method}.struct_orig,spatial_method{ispatial_method}.struct_null,'*','Color',cols(ispatial_method));
hold on
end
plot([0,max(spatial_method{ispatial_method}.struct_orig)],[0,max(spatial_method{ispatial_method}.struct_orig)]);
legend('best frame','Diff of Gaussian','Low rank approx','common temporal','spatio-temporal','45 degree line','Location','best');
xlabel('Structure in Original');
ylabel('Sructure in Null');

save('~/Nishal/Statial_null.mat','spatial_method');