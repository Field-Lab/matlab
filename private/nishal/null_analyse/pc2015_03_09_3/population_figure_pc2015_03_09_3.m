%% analyze data012,data013,data008
% Condition strings
nConditions=3;
condDuration=1272/120;
cond_str=cell(3,1);
cond_str{1}='Original';
cond_str{2}='Cell group 1 Spatial';
cond_str{3}='OFF parasol';
interestingConditions=[1,2,3];


%% make movies
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
%%
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
nullCond = 3;

NullCells1=[4083,1804,2448,5221,1068];  % OFF
NullCells2=[];   % ON


condDuration=1272/120;
GLM_fit_link= '/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_3/data008/CellType_OFFParasol/rk1/';
nConditions=3;

LL_log=zeros(length(InterestingCell_vis_id),nConditions,2);
LL_log_glm=zeros(length(InterestingCell_vis_id),2,2);


for ref_cell_number=3:length(InterestingCell_vis_id); %11
    close all
    
    cellID=InterestingCell_vis_id(ref_cell_number)
    % Plot recorded raster
    [spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_3_light(cellID,nConditions,condDuration,cond_str,neuronPath);

    LL_log(ref_cell_number,:,:) = getLL_const_var_2015_03_09_2(spkCondColl,1);
    
     interval=4;
     [resp1,resp2 ] = framework_null_expt_2015_03_09_3(cellID,GLM_fit_link,condMov{1}-0.5,interval) % give movie between -0.5 to +0.5
     glm_resp(1).spksColl=resp1;glm_resp(2).spksColl = resp2;
     LL_log_glm(ref_cell_number,:,:) = getLL_const_var_2015_03_09_2(glm_resp,2);
     
     
end

LL_diff_orig = zeros(size(LL_log,1),2);
LL_diff_orig(:,1) = LL_log(:,1,2) - LL_log(:,1,1);
LL_diff_orig(:,2) = LL_log_glm(:,1,2) - LL_log_glm(:,1,1);



LL_diff_null= zeros(size(LL_log,1),2);
LL_diff_null(:,1) = LL_log(:,nullCond,2) - LL_log(:,nullCond,1);
LL_diff_null(:,2) = LL_log_glm(:,2,2) - LL_log_glm(:,2,1);

save(sprintf('%sDiffLL.mat',GLM_fit_link),'LL_diff_orig','LL_diff_null','InterestingCell_vis_id','LL_log','LL_log_glm');

figure;plot(LL_log(:,1,2) - LL_log(:,1,1) , LL_log(:,3,2) - LL_log(:,3,1) ,'*');hold on;plot([0,0.8],[0,0.8],'g')

figure;plot(LL_diff_orig(:,2),LL_diff_orig(:,1),'*');hold on;plot(LL_diff_null(:,2),LL_diff_null(:,1),'*');hold on;plot([0,4],[0,4],'g');
xlabel('Framework');
ylabel('Data');
gmail('bhaishahster@gmail.com','Population Null update','Off Parasol done');


%%
