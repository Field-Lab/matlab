addpath(genpath('../null_analyse/'));
addpath(genpath('../null_analyse/analyse_functions'));
%startup_null_analyse_tenessee
startup_null_analyse_bertha

%%
% Condition strings
nConditions=6;
condDuration=1200/4;
cond_str=cell(nConditions,1);
cond_str{1}='Original';
cond_str{2}='Null ';

interestingConditions=[1,2,3,4,5,6,7];
% WN (0.24)
% Null (0.24)
% Null contrast controlled (0.24)
% WN (0.36)
% Null (0.36)
% WN (0.12)
% Null (0.12)
dataRuns_OFF_null = [31,32,34,38,39,40,41];
dataRuns_ON_null = [31,33,35,38,51,40,42];
movies_OFF_null =[5,6,12,9,10,13,14];
movies_ON_null = [5,8,18,9,16,13,20];
cond_dur = [20,20,20,10,10,10,10];
%% Load Movies
%  
rawMovFrames=1200/(4);
figure;
icnt=0;
% make pixel histogram
for imov=movies_OFF_null
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2016-02-17-1/Visual/null/pc2016_02_17_1_data028/%d.rawMovie',imov),rawMovFrames,1);
    subtract_movies{3}=mean(stim,1);
    subtract_movies{3}=mean(stim,1)*0+127.5;
    movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
    movie=movie/255;
    
    icnt=icnt+1;
    
    subplot(4,2,icnt);
    qq=movie;
    hist(qq(:),20)
    title(sprintf('Movie pixel histogram %d',imov));
end



% % make movies
% interval=4;
% condMov=cell(nConditions,1);
% rawMovFrames=1200/(4);
% icnt=0;
% % make pixel histogram
% for imov=[1,2,4,10,12,5,6,8]
%     [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-11-09-1/Visual/pc2015_11_09_1_fromdata009/%d.rawMovie',imov),rawMovFrames,1);
%     subtract_movies{3}=mean(stim,1);
%     subtract_movies{3}=mean(stim,1)*0+127.5;
%     movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
%     movie=movie/255;
%     
%     icnt=icnt+1;
%     qq=permute(movie,[2,3,1]);
%     ifcnt = 0;
%     condMov{icnt}=zeros(size(qq,1),size(qq,2),size(qq,3)*interval);
%     for iframe=1:size(qq,3)
%         for irepeat=1:interval
%             ifcnt=ifcnt+1;
%             condMov{icnt}(:,:,ifcnt)=qq(:,:,iframe)+0.5; % cond mov is between 0 and 1 now!
%         end
%         
%     end
%     
% end
%% contrast map 
% set up stuff for reading real responses

WN_datafile ='/Volumes/Analysis/2015-11-09-1/streamed/data009/data009'

datarun=load_data(WN_datafile)
datarun=load_params(datarun)
cellTypeId = 1;
InterestingCell_vis_id=datarun.cell_types{cellTypeId}.cell_ids; 


% make contrast map
rawMovFrames=1200/(4);
figure;
icnt=0;
cMap = cell(8,1);
h=figure('Color','w');
for imov=movies_ON_additivity
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-11-09-1/Visual/pc2015_11_09_1_fromdata009/%d.rawMovie',imov),rawMovFrames,1);
    subtract_movies{3}=mean(stim,1);
    subtract_movies{3}=mean(stim,1)*0+127.5;
    movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
    movie=movie/255;
    
    icnt=icnt+1;
    
    subplot(3,2,icnt);
    qq=movie;
    
    cMap{icnt}=contrastMap(qq);
    %if(icnt ==1 ||icnt==3)
        imagesc(repelem(cMap{icnt}(:,end:-1:1)',20,20));
        title(sprintf('cMap: %d',icnt));
        colormap gray;axis image
        caxis([0,0.5]);set(gca,'xTick',[]);set(gca,'yTick',[]);
        plot_rf_fit_nishal(datarun, InterestingCell_vis_id,'magnify',20,'fill_color',[0,1,1],'alpha',0.3,'fill',false,'edge',true,'labels',false,'edge_color',[1,0,0]);
    %end
    
%     if(icnt==2)
%   
%         h=imagesc(repelem(( cMap{icnt}(:,end:-1:1)-cMap{1}(:,end:-1:1))'./cMap{1}(:,end:-1:1)',20,20));
%         title(sprintf('cMap: %d, compared to 1',icnt));
%         caxis([-1,1]);set(gca,'xTick',[]);set(gca,'yTick',[]);
%         colormap gray;axis image
%         plot_rf_fit_nishal(datarun, InterestingCell_vis_id,'magnify',20,'fill_color',[0,1,1],'alpha',0.3,'fill',false,'edge',true,'labels',false,'edge_color',[1,0,0]);
%     end
%     
%     if(icnt==4 || icnt==5 ||icnt==6)
%         imagesc(repelem(( cMap{icnt}(:,end:-1:1)-cMap{3}(:,end:-1:1) )'./cMap{3}(:,end:-1:1)',20,20));
%         title(sprintf('cMap: %d, compared to 3',icnt));
%         caxis([-1,1]);axis image
%         colormap gray;set(gca,'xTick',[]);set(gca,'yTick',[]);
%         plot_rf_fit_nishal(datarun, InterestingCell_vis_id,'magnify',20,'fill_color',[0,1,1],'alpha',0.3,'fill',false,'edge',true,'labels',false,'edge_color',[1,0,0]);
%     end
    % caxis([3,6]);
    colorbar
    axis image
    
end

% s=hgexport('readstyle','cMap');
% hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_09_23_0/data000/cMap.eps'),s);

mask = (cMap{2}-cMap{1})~=0;
imagesc(mask);


%% additivity experiment 

dataRuns =dataRuns_ON_null;

WN_datafile = '/Volumes/Analysis/2016-02-17-1/d28_51-norefit/data028/data028';

datarun=load_data(WN_datafile)
datarun=load_params(datarun)

cellTypeId = 1
InterestingCell_vis_id=datarun.cell_types{cellTypeId}.cell_ids; 
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);

nConditions=1;

cols='krkrkrkrkrkrkkrkrkrkr';
spkCondColl=cell(length(dataRuns),1);
clear output
output(2).s=0;
for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID= InterestingCell_vis_id(ref_cell_number)
        
    for idata=1:length(dataRuns)
    Null_datafile = sprintf('/Volumes/Analysis/2016-02-17-1/d28_51-norefit/data0%02d',dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
    [spkColl,spkCondColl{idata}]=plot_raster_script_pc2016_01_05_3_light_photons(cellID,nConditions,cond_dur(idata),cond_str,neuronPath);
    end
    
%     
%     h=figure;
%     shift=0;
%     for irun = [4,5]%length(dataRuns)
%     shift =shift-max(spkCondColl{irun}.yPoints);
%     plot(spkCondColl{irun}.xPoints/20000,spkCondColl{irun}.yPoints + shift,cols(irun));
%     hold on;
%     end
%     
%     set(gca,'yTick',[]);
%     ylim([shift,0]);
%     InterestingCell_vis_id(ref_cell_number)
%     pause(1/120);


% structure in rasters
   convolve=150;
    binSz = 1/1200;len=12000;
    for conda=1:7
    realResp = makeSpikeMat(spkCondColl{conda}.spksColl, binSz,len);
    [PSTH_reca,time]=calculate_psth_fcn2(convolve,binSz,len,realResp); 
    
    cond_var = [sqrt(var(PSTH_reca(0.25*end:end)))];
    cond_mean = [sqrt(mean(PSTH_reca(0.25*end:end)))];
    output(ref_cell_number).cellID=cellID;
    output(ref_cell_number).condition(conda).PSTH = PSTH_reca(0.25*end:end);
    output(ref_cell_number).condition(conda).cond_var = cond_var;
    output(ref_cell_number).condition(conda).cond_mean = cond_mean;
    output(ref_cell_number).condition(conda).Tlen = cond_dur(conda);
    end
  
   
% correlation in structure in rasters

output(ref_cell_number).corr_direct  = zeros(7,7);
output(ref_cell_number).corr_reverse = zeros(7,7);
for conda=1:3
    for condb=1:3
        PSTH1 =  output(ref_cell_number).condition(conda).PSTH;
        PSTH2 =  output(ref_cell_number).condition(condb).PSTH;
        output(ref_cell_number).corr_direct(conda,condb) = R_2_value(PSTH1',PSTH2'); 
        output(ref_cell_number).corr_reverse(conda,condb) = R_2_value(PSTH1',PSTH2(end:-1:1)'); 
    end
end

for conda=4:5
    for condb=4:5
        PSTH1 =  output(ref_cell_number).condition(conda).PSTH;
        PSTH2 =  output(ref_cell_number).condition(condb).PSTH;
        output(ref_cell_number).corr_direct(conda,condb) = R_2_value(PSTH1',PSTH2'); 
        output(ref_cell_number).corr_reverse(conda,condb) = R_2_value(PSTH1',PSTH2(end:-1:1)'); 
    end
end

for conda=6:7
    for condb=6:7
        PSTH1 =  output(ref_cell_number).condition(conda).PSTH;
        PSTH2 =  output(ref_cell_number).condition(condb).PSTH;
        output(ref_cell_number).corr_direct(conda,condb) = R_2_value(PSTH1',PSTH2'); 
        output(ref_cell_number).corr_reverse(conda,condb) = R_2_value(PSTH1',PSTH2(end:-1:1)'); 
    end
end

end

     if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2016_02_17_1/d_null/CellType_%s',datarun.cell_types{cellTypeId}.name)))
         mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2016_02_17_1/d_null/CellType_%s',datarun.cell_types{cellTypeId}.name));
     end
    s=hgexport('readstyle','raster');
    hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2016_02_17_1/d_null/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
    save(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2016_02_17_1/d_null/CellType_%s/summary_of_responses.mat',datarun.cell_types{cellTypeId}.name),'output','WN_datafile','cellTypeId','dataRuns');

 %% 
    ON = load('/Volumes/Lab/Users/bhaishahster/analyse_2016_02_17_1/d_null/CellType_ON Parasol/summary_of_responses.mat');
    OFF = load('/Volumes/Lab/Users/bhaishahster/analyse_2016_02_17_1/d_null/CellType_OFF Parasol/summary_of_responses.mat');
    
    figure;
    data=ON;
    for icell=1:length(data.output)
    orig_var(icell) = data.output(icell).condition(1).cond_var;
    null_var(icell) = data.output(icell).condition(3).cond_var;
    end
    histogram(null_var(orig_var>0.007)./orig_var(orig_var>0.007),100);
    %scatter(orig_var,null_var);
    hold on;
    
    data=OFF;
    for icell=1:length(data.output)
    orig_var(icell) = data.output(icell).condition(1).cond_var;
    null_var(icell) = data.output(icell).condition(3).cond_var;
    end
    histogram(null_var(orig_var>0.007)./orig_var(orig_var>0.007),100);
     %scatter(orig_var,null_var);
    p