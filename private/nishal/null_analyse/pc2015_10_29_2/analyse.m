addpath(genpath('../null_analyse/'));
addpath(genpath('../null_analyse/analyse_functions'));
%startup_null_analyse_tenessee
startup_null_analyse_bertha

%%
% Condition strings
nConditions=8;
condDuration=1200/4;
cond_str=cell(nConditions,1);
cond_str{1}='Original';
cond_str{2}='Null ';

interestingConditions=[1,2,3,4,5,6,7];

dataRuns_OFF_additivity = [3,4,6,7,9,11,13];
dataRuns_ON_additivity = [3,5,6,8,10,12,13];
movies_OFF_addivitiy =[1,2,5,6,10,14,13];
movies_ON_additivity = [1,4,5,8,12,16,13];
%% Load Movies
% 
rawMovFrames=1200/(4);
figure;
icnt=0;
% make pixel histogram
for imov=[16,18,19,21,23,24,26]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-10-06-0/Visual/null/%d.rawMovie',imov),rawMovFrames,1);
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
%     [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-09-23-0/Visual/2015-09-23-0_data017/%d.rawMovie',imov),rawMovFrames,1);
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
% 
% % make contrast map
% rawMovFrames=1200/(4);
% figure;
% icnt=0;
% cMap = cell(8,1);
% h=figure('Color','w');
% for imov=[1,2,4,10,12,5,6,8]
%     [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-09-23-0/Visual/2015-09-23-0_data017/%d.rawMovie',imov),rawMovFrames,1);
%     subtract_movies{3}=mean(stim,1);
%     subtract_movies{3}=mean(stim,1)*0+127.5;
%     movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
%     movie=movie/255;
%     
%     icnt=icnt+1;
%     
%     subplot(34,2,icnt);
%     qq=movie;
%     
%     cMap{icnt}=contrastMap(qq);
%     
%     imagesc(cMap{icnt});
%     caxis([3,6]);
%     colorbar
%     axis image
%     title(sprintf('cMap: %d',imov));
% end
% 
% s=hgexport('readstyle','cMap');
% hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_09_23_0/data000/cMap.eps'),s);



%% additivity experiment -  OFF

dataRuns =dataRuns_OFF_additivity;

WN_datafile = '/Volumes/Analysis/2015-10-29-2/d00_36-norefit/data001/data001';


datarun=load_data(WN_datafile)
datarun=load_params(datarun)

cellTypeId = 2
InterestingCell_vis_id=datarun.cell_types{cellTypeId}.cell_ids; 
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);

condDuration=10;
nConditions=1;

cols='rkrkrkrkrkrkkrkrkrkr';
spkCondColl=cell(length(dataRuns),1);

for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID=InterestingCell_vis_id(ref_cell_number)
        
    for idata=1:length(dataRuns)
    Null_datafile = sprintf('/Volumes/Analysis/2015-10-29-2/d00_36-norefit/data0%02d',dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
    [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    end
    
    
    h=figure;
    for irun = 1:length(dataRuns)
    plot(spkCondColl{irun}.xPoints/20000,spkCondColl{irun}.yPoints - (irun-1)*30,cols(irun));
    hold on;
    end
    set(gca,'yTick',[]);
    ylim([-(length(dataRuns)-1)*30,30]);
    InterestingCell_vis_id(ref_cell_number)
    
    if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_2/d_add/CellType_%s',datarun.cell_types{cellTypeId}.name)))
        mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_2/d_add/CellType_%s',datarun.cell_types{cellTypeId}.name));
    end
   s=hgexport('readstyle','raster');
   hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_2/d_add/CellType_%s/CellID_%d.eps',datarun.cell_types{cellTypeId}.name,InterestingCell_vis_id(ref_cell_number)),s);
end

