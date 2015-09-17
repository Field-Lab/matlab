%% This should be STA from WN stimulus!

% Use GLM fit parameters
if(exist('null_filter','var'))
    if(strcmp(null_filter,'STA'))
    % Use WN STA
k=WNSTA(:,:,end:-1:1);
        for iframe=1:size(k,3) % Flipping? Doubt!!
            k(:,:,iframe)=k(:,:,iframe)';
        end
xcoords=1:size(WNSTA,1);
ycoords=1:size(WNSTA,2);
    end
    
    if(strcmp(null_filter,'actual'))
    k=fittedGLM.linearfilters.Stimulus.Filter(:,:,end:-1:1);
    xcoords=fittedGLM.linearfilters.Stimulus.y_coord;
    ycoords=fittedGLM.linearfilters.Stimulus.x_coord;
    end
else
    % Use WN STA
k=WNSTA(:,:,end:-1:1);
        for iframe=1:size(k,3) % Flipping? Doubt!!
            k(:,:,iframe)=k(:,:,iframe)';
        end
xcoords=1:size(WNSTA,1);
ycoords=1:size(WNSTA,2);
end



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
cell_params2.type_name_inp='cell';
% cell_params.cell_list=[3888,2825,1820,4129, 5346,5671,5161,1278, 3828,3574,4036,3572, 503,560,797,1009,487,181,901]; % if type_name_inp = 'userCellList' 
if(exist('null_filter','var'))
    if(strcmp(null_filter,'STA'))
 cell_params2.use_fits=2; % Clipped
    end
    
    if(strcmp(null_filter,'actual'))
   cell_params2.use_fits=0; % Raw
    end
else
 cell_params2.use_fits=2; % Clipped
end

if(exist('sta_spatial_method','var'))
cell_params2.sta_spatial_method=sta_spatial_method;%1,2 ,3,4
else
cell_params2.sta_spatial_method=4;
end

cell_params2.STAlen=14;
cell_params2.stas=stas_big2;

mov_params2=struct();
if(exist('mov_type_null_touse','var'))
    mov_params2.mov_type=mov_type_null_touse
else
mov_params2.mov_type='bw'
end

mov_params2.mean=0.5*255;
mov_params2.deviation=0.48*255;
%mov_params2.scaling_loss=0.01; % a number in [0,1], fraction of values that is changed by scaling.
mov_params2.stixel=10;

mov_params2.post_process_method = 'scale'; % or, 'stretch'
mov_params2.scale = 0.48/0.48;
mov_params2.interval=interval;
if(exist('interval','var'))
mov_params2.interval=interval; %'/Volumes/Analysis/stimuli/white-noise-xml/BW-10-1-0.48-11111-32x32.xml';
else
 mov_params2.interval=1;
end
mov_params2.movie_time=120*10/mov_params2.interval;

if(exist('mdf_file','var'))
mov_params2.mdf_file=mdf_file; %'/Volumes/Analysis/stimuli/white-noise-xml/BW-10-1-0.48-11111-32x32.xml';
end


if(exist('solver_use','var'))
solver=solver_use; %'/Volumes/Analysis/stimuli/white-noise-xml/BW-10-1-0.48-11111-32x32.xml';
else
solver=8;
end

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
if(exist('code_version','var'))
    if(strcmp(code_version,'old'))
        write_movie_idx(destination_mat,movie_full,mov_idx);
    else
        write_movie_idx(destination_mat,movie_full,mov_idx,mov_params2.stixel);
    end
else
    write_movie_idx(destination_mat,movie_full,mov_idx,mov_params2.stixel);
end


display(sprintf('Movie Length %d',size(movie_full,3)));
raw_mov_len = size(movie_full,3);