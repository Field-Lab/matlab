addpath(genpath('../null_analyse/'));
addpath(genpath('../null_analyse/analyse_functions'));
%startup_null_analyse_tenessee
startup_null_analyse_bertha

%%
% Condition strings
nConditions=6;
condDuration=12;
cond_str=cell(6,1);
cond_str{1}='Original';
cond_str{2}='Null for On Parasol';
cond_str{3}='Original';
cond_str{4}='Null for Off Parasol';
cond_str{5}='Original';
cond_str{6}='Null for few cells';
interestingConditions=[1,2,4,6];
%% Load Movies
rawMovFrames=8640;
[stim,height,width,header_size] = get_raw_movie('/Volumes/Analysis/nishal/pc2014_11_24_3_data012/18.rawMovie',rawMovFrames,1);
subtract_movies{3}=mean(stim,1);
subtract_movies{3}=mean(stim,1)*0+127.5;
movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
movie=movie/255;

movStart=12*120*[0:6];
condMovies=cell(6,1);
for icond=1:6
    condMovies{icond}=movie(movStart(icond)+1:movStart(icond+1),:,:);
end

figure
for icond=1:6
    subplot(3,2,icond);
    qq=condMovies{icond}(120:12*120-120,:,:);
    hist(qq(:),20)
    title(sprintf('Condition %d',icond));
end


figure('Color','w');
icnt=0;
for icondi=[2,4,6]
    icnt=icnt+1;
    subplot(3,1,icnt);
    aaa=condMovies{icondi}(240:1440-240,:,:)-condMovies{1}(240:1440-240,:,:);
    hist(aaa(:),50);
    title('Difference in pixel value compared to original movie');
end


%%


WN_datafile = 'nishal/2014-11-24-3/data012/data012';
Null_datafile = '/Volumes/Analysis/2014-11-24-3/data014';
WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/nishal/2014-11-24-3/data012';
imov=14;

datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[1]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end 
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);

for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
     [spkColl,spkCondColl]=plot_raster_script(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str);
     
    ref_cell_number
   
   plot_mosaic(datarun,InterestingCell_vis_id,ref_cell_number)
    %testsuite_prediction
  [timeLogData,psthData] = psth_variability(spkCondColl,nConditions,condMovies,cond_str,InterestingCell_vis_id,imov,ref_cell_number,interestingConditions);
    pause
end



%% Fitted v/s actual SD ? 


WN_datafile = 'nishal/2014-11-24-3/data012/data012';
Null_datafile = '/Volumes/Analysis/2014-11-24-3/data014';
WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/nishal/2014-11-24-3/data012';
imov=14;

datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[1]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end 
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);



fitsd=cell(nConditions,1);

for ref_cell_number=1:length(InterestingCell_vis_id); %11
  
    ref_cell_number
   
    %genClip=linear_output(WN_datafile,InterestingCell_vis_id,ref_cell_number,condMovies,'clip',nConditions,cond_str);
    genFit=linear_output(WN_datafile,InterestingCell_vis_id,ref_cell_number,condMovies,'fit',nConditions,cond_str);
    
    %clipsd=[clipsd;sqrt(var(genClip))];
    for icond=1:nConditions
    fitsd{icond}=[fitsd{icond};sqrt(var(genFit{icond}))];
    end
end

clipsd=cell(nConditions,1);
for ref_cell_number=1:length(InterestingCell_vis_id);
    ref_cell_number
    genClip=linear_output(WN_datafile,InterestingCell_vis_id,ref_cell_number,condMovies,'clip',nConditions,cond_str);
    for icond=1:nConditions
    clipsd{icond}=[clipsd{icond};sqrt(var(genClip{icond}))];
    end
end

%% Plots
figure('Color','w')

for icond=2:4
    subplot(2,2,icond);
    
    plot(fitsd{1},fitsd{icond},'.r');
    hold on
    plot(clipsd{1},clipsd{icond},'.b');
    xlim([0,0.7]);
    ylim([0,0.7]);
    hold on
    plot(0:0.1:0.7,0:0.1:0.7,'r');
    xlabel('Condition 1 sd');
    ylabel(sprintf('Condition %s sd',cond_str{icond}));
    
    legend('Fitted STA','Clipped STA','y=x line','Location','best');
    title(sprintf('On Parasol Cond1 v/s %s',cond_str{icond}));
end

%%



WN_datafile = 'nishal/2014-11-24-3/data012/data012';
Null_datafile = '/Volumes/Analysis/2014-11-24-3/data014';
WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/nishal/2014-11-24-3/data012';
imov=14;

datarun=load_data(WN_datafile)
datarun=load_params(datarun)


cellTypeId=[2]; % 1 for On Parasols, 2 for Off parasols
InterestingCell_vis_id=[];
for icellType=cellTypeId
    icellType 
    InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
end 
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);


clipsd2=cell(nConditions,1);
for ref_cell_number=1:length(InterestingCell_vis_id);
    ref_cell_number
    genClip=linear_output(WN_datafile,InterestingCell_vis_id,ref_cell_number,condMovies,'clip',nConditions,cond_str);
    for icond=1:nConditions
    clipsd2{icond}=[clipsd2{icond};sqrt(var(genClip{icond}))];
    end
end


fitsd2=cell(nConditions,1);

for ref_cell_number=1:length(InterestingCell_vis_id); %11
  
    ref_cell_number
   close all
  
    try
    %genClip=linear_output(WN_datafile,InterestingCell_vis_id,ref_cell_number,condMovies,'clip',nConditions,cond_str);
    genFit=linear_output(WN_datafile,InterestingCell_vis_id,ref_cell_number,condMovies,'fit',nConditions,cond_str);
    
    %clipsd=[clipsd;sqrt(var(genClip))];
    for icond=1:nConditions
    fitsd2{icond}=[fitsd2{icond};sqrt(var(genFit{icond}))];
    end
    catch
    for icond=1:nConditions
    fitsd2{icond}=[fitsd2{icond};-1];
    end
        
    end
end

%% Plots
figure('Color','w')

for icond=2:6
    subplot(3,2,icond);
    
    plot(fitsd2{1},fitsd2{icond},'.r');
    hold on
    plot(clipsd2{1},clipsd2{icond},'.b');
    xlim([0,0.7]);
    ylim([0,0.7]);
    hold on
    plot(0:0.1:0.7,0:0.1:0.7,'r');
    xlabel('Condition 1 sd');
    ylabel(sprintf('Condition %s sd',cond_str{icond}));
    
    legend('Fitted STA','Clipped STA','y=x line','Location','best');
    title(sprintf('Off Parasol Cond1 v/s %s',cond_str{icond}));
end

save('/Volumes/Analysis/nishal/Null/pc2014_11_24_3.mat');