%% This should be STA from WN stimulus!

% Use GLM fit parameters
% k=fittedGLM.linearfilters.Stimulus.Filter;
% xcoords=fittedGLM.linearfilters.Stimulus.y_coord;
% ycoords=fittedGLM.linearfilters.Stimulus.x_coord;

% Use WN STA
k=WNSTA(:,:,end:-1:1);
        for iframe=1:size(k,3) % Flipping? Doubt!!
            k(:,:,iframe)=k(:,:,iframe)';
        end
xcoords=1:size(WNSTA,1);
ycoords=1:size(WNSTA,2);

datafile = 'load_from_cell_params';
stas_big{1}=zeros(32,32,3,30); %zeros(32,64,3,30);
stas_big{1}(xcoords,ycoords,1,1:end)=k;
stas_big{1}(xcoords,ycoords,2,1:end)=k;
stas_big{1}(xcoords,ycoords,3,1:end)=k;
stas_big2{1}=stas_big{1};%+rand(size(stas_big{1}))*0.06;

figure;
for itime=1:30
    itime
imagesc(sum(stas_big2{1}(:,:,:,itime),3))
colormap gray
colorbar
caxis([3*min(min(min(min(stas_big2{1}(:,:,:,:))))),3*max(max(max(max(stas_big2{1}(:,:,:,:)))))]);
pause(1/120); 
end

%% Replace 'stas_big2' with something WN-STA  
destination_mat='/Volumes/Lab/Users/bhaishahster/Spatial_null/Figures_EJ';

cell_params2=struct();
cell_params2.type_name_inp='nc2';%'userCellList';
cell_params2.cell_list=[]%[3888,2825,1820,4129, 5346,5671,5161,1278, 3828,3574,4036,3572, 503,560,797,1009,487,181,901]; % if type_name_inp = 'userCellList' 
cell_params2.STAlen=14;
cell_params2.sta_spatial=[];%sprintf('%s/stas_spatial_test_suite.mat',destination_mat);
cell_params2.use_fits=use_fit_var; % 2, 0,0,2
cell_params2.sta_spatial_method=sta_spatial_method_var;%1,2 ,3,4
cell_params2.stas=stas_big2;
% cell_params.sta_spatial_method = 1 for just using 4th frame, 2 is for fitting spatial STA. 
% Use cell_params.use_fits=2 (clipping) if cell_params.sta_spatial_method = 1 and 
% use cell_params.use_fits=0 (no processing of STA) if
% cell_params.sta_spatial_method = 2;

mov_params2=struct();
mov_params2.mov_type='bw';
mov_params2.movie_time=120*10;
mov_params2.mean=0.5*255;
mov_params2.deviation=0.48*255;
mov_params2.scaling_loss=0.01; % a number in [0,1], fraction of values that is changed by scaling.


solver=4; % Solver 4 used for spatial nulling!
[mov_orignial,mov_modify_new]=null_space_movie2(datafile,cell_params2,mov_params2,solver);

% Write movies
movies{1}=mov_orignial;
movies{2}=mov_modify_new;


movie_list= [1,2];
movie_full=zeros(size(movies{1},1),size(movies{1},2),size(movies{1},3)*length(movie_list));
icnt=1;
for imov=movie_list
movie_full(:,:,(icnt-1)*size(movies{imov},3)+1:icnt*size(movies{imov},3))=movies{imov};
icnt=icnt+1;
end
mov_idx=18;
write_movie_idx(destination_mat,movie_full,mov_idx,10);
display(sprintf('Movie Length %d',size(movie_full,3)));

%%