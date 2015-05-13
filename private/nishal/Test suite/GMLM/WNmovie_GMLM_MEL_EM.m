%% GMLM with MEL and EM for normal cells at coarser resolution .. 

location_of_git_repo='/home/vision/Nishal/matlab/';
addpath(genpath([location_of_git_repo '/private/nora']));

% Java library

addpath(('~/Nishal/matlab/private/nishal/create_act2'));
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');
addpath(genpath('~/Nishal/matlab/private/nishal/create_act_2/'));
addpath(genpath('~/Nishal/matlab/private/nishal/Test suite/GLM/'));
addpath(genpath('~/Nishal/matlab/private/nishal/Test suite/GLM2/'));
addpath(genpath('~/Nishal/matlab/private/nishal'));
addpath(genpath('~/Nishal/matlab/code'));
%% Dataset details
WN_datafile = '2015-03-09-2/streamed/data038/data038';
WN_datafile_short='2015-03-09-2/streamed/data038/data038';
movie_xml = 'BW-8-2-0.48-11111-40x40';
stim_length=1800;% in seconds
cellID = 2042;
cell_glm_fit = '/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data038/CellType_OFF parasol/CellID_2042.mat';

%% Load stimuli


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);


%
stim_description=movie_xml;
xml_file=['/Volumes/Analysis/stimuli/white-noise-xml/' stim_description '.xml'];
dashes=find(stim_description=='-');
StimulusPars.type=stim_description(1:dashes(1)-1);
StimulusPars.pixelsize = str2double(stim_description(dashes(1)+1:dashes(2)-1));
StimulusPars.refreshrate = str2double(stim_description(dashes(2)+1:dashes(3)-1));
StimulusPars.RNG = str2double(stim_description(dashes(3)+1:dashes(4)-1));
% try
%     StimulusPars.height = str2double(stim_description(dashes(5)+1:dashes(5)+2)); 
%     StimulusPars.width = str2double(stim_description(dashes(5)+4:dashes(5)+5));
% catch
%     StimulusPars.height = 32; 
%     StimulusPars.width = 64;
% end
StimulusPars.tstim = 1/120;
fitframes=stim_length*120; % seconds * 120 frames per second / interval


disp('Loading Stimulus Movies')
[temp_fitmovie,height,width,~,~] = get_movie(xml_file, datarun.triggers, fitframes/StimulusPars.refreshrate);
temp_fitmovie=permute(temp_fitmovie,[2 1 3 4]);
fitmovie_color=zeros(width,height,3,fitframes);

try
    StimulusPars.height = str2double(stim_description(dashes(5)+1:dashes(5)+2)); 
    StimulusPars.width = str2double(stim_description(dashes(5)+4:dashes(5)+5));
catch
    StimulusPars.height =height; 
    StimulusPars.width = width;
end


for i=1:fitframes
    fitmovie_color(:,:,:,i)=temp_fitmovie(:,:,:,ceil(i/StimulusPars.refreshrate));
end

clear temp_fitmovie i 

%% Load spikes 
    master_idx         = find(datarun.cell_ids == cellID);
    
      
        % Spike loading
        spikes=datarun.spikes{master_idx};
    
        % make STA 3D ? 
        %glm_cellinfo.WN_STA = squeeze(sum(glm_cellinfo.WN_STA,3)); % Doubt!!!!!!!
        clear cell_savename
        
        % Align the spikes and the movies;
        spikes_adj=spikes;
        n_block=0;
        for i=1:(length(datarun.triggers)-1)
            actual_t_start=datarun.triggers(i);
            supposed_t_start=n_block*100/120;
            idx1=spikes > actual_t_start;
            idx2=spikes < datarun.triggers(i+1);
            spikes_adj(find(idx2.*idx1))=spikes(find(idx2.*idx1))+supposed_t_start-actual_t_start;
            n_block=n_block+1;
        end
        clear spikes
        spike.home=spikes_adj;
        clear spikes_adj;
        %%
        spksGen = zeros(stim_length*120,1);
        for ispike=1:length(spike.home)
            spksGen(floor(spike.home(ispike)*120)+1)=1;
        end
        spksGen = spksGen(1:stim_length*120);
        
        Filtlen=30;
        
         %% find significant mask - 1
        idx=1:length(datarun.cell_ids);
        matlab_cell_ids = idx(datarun.cell_ids==cellID);
        stas=datarun.stas.stas(matlab_cell_ids);
        
        % Load STAs
        
        stas_new=cell(length(stas),1);
        for icell=1:length(stas)
            st_temp=zeros(size(stas{1},2),size(stas{1},1),1,size(stas{1},4)); % DOUBT .. Could be a reason for things to fail!!!!!
            for itime=1:size(stas{1},4)
                st_temp(:,:,:,itime)=mean(stas{icell}(:,:,:,end-itime+1),3)'; % DOUBT .. Could be a reason for things to fail!!!!!
            end
            stas_new{icell}=st_temp;
        end
        stas=stas_new;
        
        % Used in movie post process
        cell_params.STAlen=14;
        [stas_clipped,totalMaskAccept2,CellMasks]= clipSTAs(stas,cell_params);
        
        mov = fitmovie_color-0.5;
        iisp = 1:size(mov,4);
        iisp = iisp(spksGen~=0 & iisp' > 31);
        STA_recalc = zeros(40,40,3,30);
        for ilen=1:length(iisp)
        STA_recalc = STA_recalc + mov(:,:,:,iisp(ilen)-29:iisp(ilen));
        end
        STA_recalc =STA_recalc/numel(iisp);
        
        figure; 
        subplot(1,2,1);
        imagesc(mean(STA_recalc(:,:,:,24),3));
        colormap gray
        subplot(1,2,2);
        imagesc(totalMaskAccept2);
        colormap gray
        %tf ? 
        

      tf=squeeze(mean(mean(mean(STA_recalc.*repmat(totalMaskAccept2,[1,1,3,30]),3),1),2));
      %tf=tf/max(abs(tf));
      %idx=1:length(tf);
      tf=tf(end:-1:1);
      %tf=tf.*double(idx<15)';
     
      figure;
        plot(tf);
     
 mov=squeeze(mean(mov,3));
 maskedMovdd= filterMov(mov,totalMaskAccept2,squeeze(tf));
 maskedMov2dd=[maskedMovdd;ones(1,size(maskedMovdd,2))];
% [fitGMLM,output] = fitGMLM_afterSTC_simplified(binnedResponses,maskedMov,7,4);
 
% x_coord=[26-13:26+13]; y_coord = [26-13:26+13];
%% Load old GLM fits
% fitGLM = load (cell_glm_fit);
% fitGLM=fitGLM.fittedGLM;
% tf = fitGLM.linearfilters.Stimulus.time_rk1;
% figure;
% plot(tf);
% x_coord = fitGLM.linearfilters.Stimulus.x_coord;
% y_coord = fitGLM.linearfilters.Stimulus.y_coord;
% totalMaskAccept = zeros(40,40);
% totalMaskAccept(x_coord,y_coord)=1;
% 
%  
%         mov = fitmovie_color-0.5;
%         iisp = 1:size(mov,4);
%         iisp = iisp(spksGen~=0 & iisp' > 31);
%         STA_recalc = zeros(40,40,3,30);
%         for ilen=1:length(iisp)
%         STA_recalc = STA_recalc + mov(:,:,:,iisp(ilen)-29:iisp(ilen));
%         end
%         STA_recalc =STA_recalc/numel(iisp);
%         
%         figure; 
%         subplot(1,2,1);
%         imagesc(mean(STA_recalc(:,:,:,26),3));
%         colormap gray
%         subplot(1,2,2);
%         imagesc(totalMaskAccept);
%         colormap gray
%         
%         figure;
%         imagesc(mean(STA_recalc(x_coord,y_coord,:,26),3));colormap gray
%         
%  mov=squeeze(mean(mov,3));
%  maskedMovdd= filterMov(mov,totalMaskAccept,squeeze(tf));
%  maskedMov2dd=[maskedMovdd;ones(1,size(maskedMovdd,2))];

 %% EM like Max Expected Likelihood .. 
 interval=1;
%  idx = [1:end]
 binnedResponsesbigd = spksGen;
 mov_use=maskedMov2dd;
 filteredStimDim=size(mov_use,1);
 nSU = 3;
 [fitGMLM,output] = fitGMLM_MEL_EM(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval);  
 %%
 [fitGMLM_full2,output]= fitGMLM_full(fitGMLM,spike.home,mov_use);
figure;
 plot(fitGMLM_full2.hist.hexpanded)
  %% Show learned filters;
  mask = totalMaskAccept2;
sta_dim1 = size(mask,1);
sta_dim2 = size(mask,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask));

u_spatial_log = zeros(40,40,nSU);

figure;
for ifilt=1:nSU
subplot(2,2,ifilt)
u_spatial = reshape_vector(fitGMLM_full2.Linear.filter{ifilt}(1:length(masked_frame)),masked_frame,indexedframe);
imagesc(u_spatial(x_coord,y_coord));
colormap gray
colorbar
title(sprintf('GMLM Filter: %d',ifilt));
axis square
u_spatial = u_spatial.*totalMaskAccept;
u_spatial_log(:,:,ifilt) = u_spatial;
end

[v,I] = max(u_spatial_log,[],3);

xx=I.*(v>0.2);
xx=xx(x_coord,y_coord);
figure;
imagesc(xx);



%% Response prediction 

%% Load different movies
nConditions=6;
% make movies
interval=2;
condMov=cell(nConditions,1);
rawMovFrames=1272/(2);
icnt=0;
% make pixel histogram
for imov=[1,2,4,6,8,10]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-03-09-2/Visual/pc2015_03_09_2_data038/%d.rawMovie',imov),rawMovFrames,1);
    subtract_movies{3}=mean(stim,1);
    subtract_movies{3}=mean(stim,1)*0+127.5;
    movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
    movie=movie/255;
    
    icnt=icnt+1;
    qq=permute(movie,[2,3,1]);
    ifcnt = 0;
    condMov{icnt}=zeros(size(qq,1),size(qq,2),size(qq,3)*interval);
    for iframe=1:size(qq,3)
        for irepeat=1:interval
            ifcnt=ifcnt+1;
            condMov{icnt}(:,:,ifcnt)=double(qq(:,:,iframe)); % cond mov is between 0 and 1 now!
        end
        
    end
    
end

%% Do predictions - full
% 
% % Load recorded response
% Null_datafile = '/Volumes/Analysis/2015-03-09-2/data042-from-data038_streamed_nps';
% neuronPath = [Null_datafile,sprintf('/data042-from-data038_streamed_nps.neurons')];
% condDuration=10.6;
% nConditions=6;
% cond_str=[];
% [spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_2_light(cellID,nConditions,condDuration,cond_str,neuronPath);
% 
% 
% % Make predictions
% pred=cell(nConditions,1);
% for  icond=1:nConditions
% movd = condMov{icond};
%  maskedMov= filterMov(movd,mask,squeeze(tf));
%  maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
% nTrials=50;
%  pred{icond}= predictGMLM_full(fitGMLM_full2,maskedMov2,nTrials)';
% end
% 
% plot_record_prediction(spkCondColl,pred)

%% Do predictions - no feedback

% Load recorded response
Null_datafile = '/Volumes/Analysis/2015-03-09-2/data042-from-data038_streamed_nps';
neuronPath = [Null_datafile,sprintf('/data042-from-data038_streamed_nps.neurons')];
condDuration=10.6;
nConditions=6;
cond_str=[];
[spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_2_light(cellID,nConditions,condDuration,cond_str,neuronPath);


% Make predictions
pred=cell(nConditions,1);
for  icond=1:nConditions
movd = condMov{icond};
 maskedMov= filterMov(movd,mask,squeeze(tf));
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
nTrials=50;
 pred{icond}= predictGMLM(fitGMLM,maskedMov2,nTrials)';
end

plot_record_prediction(spkCondColl,pred)