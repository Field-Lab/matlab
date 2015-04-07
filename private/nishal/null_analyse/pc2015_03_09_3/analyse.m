addpath(genpath('../null_analyse/'));
addpath(genpath('../null_analyse/analyse_functions'));
%startup_null_analyse_tenessee
startup_null_analyse_bertha

%%
% Condition strings
nConditions=5;
condDuration=1270/120;
cond_str=cell(6,1);
cond_str{1}='Original';
cond_str{2}='Cell group 1 Spatial';
cond_str{3}='Cell group 1 Spatio-temporal';
cond_str{4}='Cell group 2 Spatial';
cond_str{5}='On Parasol';
interestingConditions=[1,2,3,4,5];

%% Load Movies
rawMovFrames=1270/(1);
figure;
icnt=0;
for imov=[1,2,4,6,8]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-03-09-3/Visual/Null Movies/pc2015_03_09_3_data000/%d.rawMovie',imov),rawMovFrames,1);
    subtract_movies{3}=mean(stim,1);
    subtract_movies{3}=mean(stim,1)*0+127.5;
    movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
    movie=movie/255;

    icnt=icnt+1;
    
    subplot(3,2,icnt);
    qq=movie;
    hist(qq(:),20)
    title(sprintf('Movie pixel histogram %d',imov));
end

rawMovFrames=1270/(1);
figure;
icnt=0;
cMap = cell(6,1);
h=figure('Color','w');
for imov=[1,2,4,6,8]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-03-09-3/Visual/Null Movies/pc2015_03_09_3_data000/%d.rawMovie',imov),rawMovFrames,1);
    subtract_movies{3}=mean(stim,1);
    subtract_movies{3}=mean(stim,1)*0+127.5;
    movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
    movie=movie/255;

    icnt=icnt+1;
    
    subplot(3,2,icnt);
    qq=movie;
   
    cMap{icnt}=contrastMap(qq);
   
    imagesc(cMap{icnt});
    %caxis([3,6]);
    colorbar
    axis image
    title(sprintf('cMap: %d',imov));
end

   s=hgexport('readstyle','cMap');
   hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data000/cMap.eps'),s);
  
%% data003 from data000


WN_datafile = '2015-03-09-3/streamed/data000/data000';
Null_datafile = '/Volumes/Analysis/2015-03-09-3/data003-04-from-data000_streamed_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2015-03-09-3/streamed/data000/data000';
neuronPath = [Null_datafile,sprintf('/data003-04-from-data000_streamed_nps.neurons')];


datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[1]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end 

%InterestingCell_vis_id = [5568,1726,5252,3061]; % OFF 
% InterestingCell_vis_id = [6106,1835,7730]; % ON
%InterestingCell_vis_id = [4501]; % SBC

InterestingCell_vis_id=[7368,872,4114,2448,5581];
NullCells1=[7368,872,4114,2448,5581];  % OFF
NullCells2=[4863,2615,3783,872,5735];   % ON


condDuration=1270/120;
nConditions=5;
 PSTH_var_log=zeros(length(InterestingCell_vis_id),3);
for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number);
    % [spkColl,spkCondColl,h]=plot_raster_script_pc2015_02_24_2(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str,neuronPath);
 [spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_3_light(cellID,nConditions,condDuration,cond_str,neuronPath);
 
   for icond=1:3
                    rec_rast= makeSpikeMat(spkCondColl(icond).spksColl,1/120,rawMovFrames);
                    [PSTH_rec,time]=calculate_psth_fcn2(100,1/120,rawMovFrames,rec_rast);
                    PSTH_rec=PSTH_rec/norm(PSTH_rec);
                    PSTHLen=length(PSTH_rec); idx=1:PSTHLen;
                    PSTH_rec=PSTH_rec(idx>0.1*PSTHLen & idx<0.9*PSTHLen);
                    PSTH_var_log(ref_cell_number,icond)=sqrt(var(PSTH_rec));
    end 
 
 plot_mosaic_pc2015_03_09_3(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2);
   InterestingCell_vis_id(ref_cell_number)
    if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data003/CellType_%s',datarun.cell_types{cellTypeId}.name)))
    mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data003/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','ras_mos4');
   hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data003/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
  
    %testsuite_prediction
  %[timeLogData,psthData] = psth_variability(spkCondColl,nConditions,condMovies,cond_str,InterestingCell_vis_id,imov,ref_cell_number,interestingConditions);
   % pause
end



%% analyze data012,data013,data008
% Condition strings
nConditions=3;
condDuration=1272/120;
cond_str=cell(3,1);
cond_str{1}='Original';
cond_str{2}='Cell group 1 Spatial';
cond_str{3}='OFF parasol';
interestingConditions=[1,2,3];

%% data012 movies
% make movies
interval=4;
condMov=cell(nConditions,1);
rawMovFrames=1272/(interval);
icnt=0;
% make pixel histogram
for imov=[1,2,4]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-03-09-3/Visual/Null Movies/pc2015_03_09_3_data008/%d.rawMovie',imov),rawMovFrames,1);
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
            condMov{icnt}(:,:,ifcnt)=qq(:,:,iframe)+0.5; % cond mov is between 0 and 1 now!
        end
        
    end
    
end
%% data012-from-data008
WN_datafile = '2015-03-09-3/streamed/data008/data008';
Null_datafile = '/Volumes/Analysis/2015-03-09-3/data012-13-from-data008_streamed_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2015-03-09-3/streamed/data008/data008';
neuronPath = [Null_datafile,sprintf('/data012-13-from-data008_streamed_nps.neurons')];


datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[2]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end 

InterestingCell_vis_id = 1068;%[4083,1804,2448,5221,1068]; % OFF 


NullCells1=[4083,1804,2448,5221,1068];  % OFF
NullCells2=[];   % ON


condDuration=1272/120;
GLM_fit_link= '/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data008/CellType_OFFParasol/';
nConditions=3;
for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number);
    % [spkColl,spkCondColl,h]=plot_raster_script_pc2015_02_24_2(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str,neuronPath);
 [spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_3_light(cellID,nConditions,condDuration,cond_str,neuronPath);
 
 plot_mosaic_pc2015_03_09_3(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2);
   InterestingCell_vis_id(ref_cell_number)
    if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data012/CellType_%s/CellID_%d',datarun.cell_types{cellTypeId}.name,cellID)))
    mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data012/CellType_%s/CellID_%d',datarun.cell_types{cellTypeId}.name,cellID));
    end
   s=hgexport('readstyle','ras_mos4');
   hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data012/CellType_%s/CellID_%d/recorded.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
  
   
   
    % GLM predictions
    figure;
    rank=2;
    [spkCondCollGLM,h2]=plot_GLM_prediction_pc2015_03_09_3(cellID,condMov,GLM_fit_link,rank);
    plot_mosaic_pc2015_03_09_3(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2);
    
    h3=plot_record_prediction_pc2015_03_09_2(spkCondColl,spkCondCollGLM);
    plot_mosaic_pc2015_03_09_2(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2);
    s=hgexport('readstyle','ras_mos4');
    hgexport(h3,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data012/CellType_%s/CellID_%d/GLM_pred_rank%d_recorded.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number),rank),s);
  
    %testsuite_prediction
  %[timeLogData,psthData] = psth_variability(spkCondColl,nConditions,condMovies,cond_str,InterestingCell_vis_id,imov,ref_cell_number,interestingConditions);
   % pause
end

%% 
    correlation=zeros(2,6,5,3);
    nullscale_list = [1,5,10,25,50,100,2];
    origscale_list = [1,10];
    convolve_list = [50,100,150,200,400];
    cutofflow=0.2;
    cutoffhigh=0.1;
    for iorigscale=[1]%1:2
        orig_scale=origscale_list(iorigscale);
        for inullscale=[1]%1:6
            null_scale=nullscale_list(inullscale);
            figure;
            mov_scales=[orig_scale,null_scale,null_scale,1,1,1];
            [spkCondCollGLM,h7]=plot_GLM_prediction_movie_scaled_pc2015_03_09_3(cellID,condMov,GLM_fit_link,mov_scales,rank);
            
             h8=plot_record_prediction_pc2015_03_09_2(spkCondColl,spkCondCollGLM);
             plot_mosaic_pc2015_03_09_2(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2);
              s=hgexport('readstyle','ras_mos4');
             hgexport(h8,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data012/CellType_%s/CellID_%d/contrast_raster_scale_orig_%d_null_%d_rank_%d.eps',datarun.cell_types{cellTypeId}.name,cellID,orig_scale,null_scale,rank),s);
    
            h=  figure('Color','w')
            for icond=1:3
                for iconvolve=2%1:5
                   
                    convolve=convolve_list(iconvolve);
                    rec_rast=spkCondCollGLM{icond};
                    [PSTH_rec,time]=calculate_psth_fcn(convolve,fittedGLM,rec_rast);
                    PSTH_pred = PSTH_rec;
                    PSTH_pred=PSTH_pred/norm(PSTH_pred);
                    
                    rec_rast= makeSpikeMat(spkCondColl(icond).spksColl,1/1200,size(spkCondCollGLM{icond},2));
                    [PSTH_rec,time]=calculate_psth_fcn(convolve,fittedGLM,rec_rast);
                    PSTH_rec=PSTH_rec/norm(PSTH_rec);
                     
                    times2=time(1:length(time));
                    time_len = max(times2);
                    cond_times=times2>cutofflow*time_len & times2<(1-cutoffhigh)*time_len;
   
                    PSTH_pred = PSTH_pred(cond_times);
                    PSTH_rec = PSTH_rec(cond_times);
                    
                  subplot(3,1,icond);
                    plot(times2(cond_times),PSTH_pred,'k');
                    hold on;
                    plot(times2(cond_times),PSTH_rec,'r');
                    h_legend=legend('Predicted','Recorded');
                    set(h_legend,'FontSize',5,'Location','best');
                    
                    correlation(iorigscale,inullscale,iconvolve,icond) = (PSTH_pred-mean(PSTH_pred))*(PSTH_rec-mean(PSTH_rec))';
                    title(sprintf('Scale : %d Correlation %f',mov_scales(icond),correlation(iorigscale,inullscale,iconvolve,icond)));
                    [iorigscale,inullscale,icond,iconvolve,correlation(iorigscale,inullscale,iconvolve,icond)]
                    

%                     subplot(3,1,icond);
%                     plot(xcorr((PSTH_pred-mean(PSTH_pred)),(PSTH_rec-mean(PSTH_rec))),'r');
%                     h_legend=legend('Cross-Correlation');
%                     set(h_legend,'FontSize',5,'Location','best');
%                     
%                     correlation(iorigscale,inullscale,iconvolve,icond) = (PSTH_pred-mean(PSTH_pred))*(PSTH_rec-mean(PSTH_rec))';
%                     title(sprintf('Scale : %d Correlation %f',mov_scales(icond),correlation(iorigscale,inullscale,iconvolve,icond)));
%                     [iorigscale,inullscale,icond,iconvolve,correlation(iorigscale,inullscale,iconvolve,icond)]
%                     
                end
            end
             
    s=hgexport('readstyle','psth_conditions');
    hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data012/CellType_OFF Parasol/CellID_%d/contrast_psth__xcorr_scale_orig_%d_null_%d_rank_%d',cellID,orig_scale,null_scale,rank),s);
    
        end
    end
   
    % Experiments
%     figure;
%     plot(squeeze(correlation(2,:,2,:)))
%     legend('Condition 1','Condition 2','Condition 3');


%% data013
%% analyze data013,data008
% Condition strings
nConditions=2;
cond_str=cell(3,1);
cond_str{1}='Original';
cond_str{2}='Cell group 1 Spatial';
cond_str{3}='OFF parasol';
interestingConditions=[1,2,3];

%% data013 movies
% make movies
interval=4;
condMov=cell(nConditions,1);
rawMovFrames=10872/(interval);
icnt=0;
% make pixel histogram
for imov=[11,12]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-03-09-3/Visual/Null Movies/pc2015_03_09_3_data008/%d.rawMovie',imov),rawMovFrames,1);
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
            condMov{icnt}(:,:,ifcnt)=qq(:,:,iframe)+0.5; % cond mov is between 0 and 1 now!
        end
        
    end
    
end


%% data013-from-data008
WN_datafile = '2015-03-09-3/streamed/data008/data008';
Null_datafile = '/Volumes/Analysis/2015-03-09-3/data013-from-data008_s_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2015-03-09-3/streamed/data008/data008';
neuronPath = [Null_datafile,sprintf('/data013-from-data008_s_nps.neurons')];


datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[2]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end 

InterestingCell_vis_id =[4083,2448,5221,1068,1804]%[4083,1804,2448,5221,1068]; % OFF 


NullCells1=[4083,1804,2448,5221,1068];  % OFF
NullCells2=[];   % ON


condDuration=10872/120;

nConditions=2;
for ref_cell_number=2:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number);
    % [spkColl,spkCondColl,h]=plot_raster_script_pc2015_02_24_2(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str,neuronPath);
 [spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_3_light(cellID,nConditions,condDuration,cond_str,neuronPath);
 
 plot_mosaic_pc2015_03_09_3(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2);
   InterestingCell_vis_id(ref_cell_number)
    if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data013/CellType_%s/CellID_%d',datarun.cell_types{cellTypeId}.name,cellID)))
    mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data013/CellType_%s/CellID_%d',datarun.cell_types{cellTypeId}.name,cellID));
    end
   s=hgexport('readstyle','ras_mos4');
   hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data013/CellType_%s/CellID_%d/recorded.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);

    % GLM predictions
    figure;
    GLM_fit_link= '/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data008/CellType_OFFParasol';
    rank=1;
    [spkCondCollGLM,h2]=plot_GLM_prediction_pc2015_03_09_3(cellID,condMov,GLM_fit_link,rank);
    plot_mosaic_pc2015_03_09_3(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2);
    
    h3=plot_record_prediction_pc2015_03_09_2(spkCondColl,spkCondCollGLM);
    plot_mosaic_pc2015_03_09_2(datarun,InterestingCell_vis_id,ref_cell_number,NullCells1,NullCells2);
    s=hgexport('readstyle','ras_mos4');
    hgexport(h3,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data013/CellType_%s/CellID_%d/GLM_pred_rank%d_recorded.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number),rank),s);
  
      
   % find STA and STC from classication run.
%    WN_datafile = '2015-03-09-3/streamed/data008/data008';
movie_xml = 'RGB-8-4-0.48-11111';
stim_len=900;% in seconds
datarun=load_data(WN_datafile)
datarun=load_params(datarun)
%cellID=3121;
dsave=sprintf('/Volumes/Lab/Users/bhaishahster/STC_computation');
Mask=zeros(80,40);
[WNSTA,dimensions,Mask]=STA_STC_from_WNrun({cellID}, WN_datafile, movie_xml, stim_len,dsave,Mask)
InterestingCell_vis_id(ref_cell_number)
dimensions_in=dimensions;



    % null spike triggered analysis
     for condAnal=1:2 
   % Mask=zeros(80,40);
    Mask = dimensions_in.Mask;
     %dimensions_in=[];
      xx=datarun.cell_types{cellTypeId}.name;
      xx(xx==' ') = '';
      
      GLM_fit_link= '/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data008/CellType_OFFParasol/rk1';
      d_save = sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data013/CellType_%s/cond_%d/fixedSP_rk1_linear',xx,condAnal);
      [reSTA,h,h2,dimensions_out]= spike_triggred_analysis_pc2015_03_09_3(condMov,1/120,neuronPath,cellID,nConditions,condDuration,condAnal,WN_datafile,d_save,GLM_fit_link,Mask,dimensions_in);
      sa=hgexport('readstyle','reSTA');
      hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data013/CellType_%s/CellID_%d/nullSTA_mov_%d_GLMrank%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number),condAnal,rank),sa);
      
      sa=hgexport('readstyle','projections');
      hgexport(h2,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data013/CellType_%s/CellID_%d/projections_mov_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number),condAnal),sa);
     end
     
%          % null spike triggered analysis
%          dimensions_in=[];
%            Mask=zeros(80,40);
% for condAnal=1:2 
%   
%      
%       xx=datarun.cell_types{cellTypeId}.name;
%       xx(xx==' ') = '';
%       
%       GLM_fit_link= '/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data008/CellType_OFFParasol/rk1';
%       d_save = sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data013/CellType_%s/cond_%d/fixedSP_rk1_linear',xx,condAnal);
%       [reSTA,h,h2,dimensions]= spike_triggred_analysis_pc2015_03_09_3(condMov,1/120,neuronPath,cellID,nConditions,condDuration,condAnal,WN_datafile,d_save,GLM_fit_link,Mask,dimensions_in);
% 
% if(condAnal==2)
%       sa=hgexport('readstyle','reSTA');
%       hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data013/CellType_%s/CellID_%d/nullSTA_mov_%d_Mask_1_GLMrank%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number),condAnal,rank),sa);
%       
%       sa=hgexport('readstyle','projections');
%       hgexport(h2,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data013/CellType_%s/CellID_%d/projections_mov_%d_dimensions_1.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number),condAnal),sa);
% end
%      dimensions_in=dimensions;
%      Mask=dimensions.Mask;
% end

     
     
    %testsuite_prediction
  %[timeLogData,psthData] = psth_variability(spkCondColl,nConditions,condMovies,cond_str,InterestingCell_vis_id,imov,ref_cell_number,interestingConditions);
   % pause
end
