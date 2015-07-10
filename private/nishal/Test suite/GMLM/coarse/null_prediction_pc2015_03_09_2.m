%% prediction to null stimulus for GLM fits. 
% Dataset: 2015-03-09-2/data038 (WN) and data042 (null run), NDF0

%% Load dataset information

WN_datafile = '2015-03-09-2/streamed/data038/data038';
WN_datafile_short='2015-03-09-2/streamed/data038/data038';
movie_xml = 'BW-8-2-0.48-11111-40x40';
wrong_xml = 'BW-8-2-0.48-22222-40x40';
stim_length=1800;% 


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);



%% Load different null movies
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

%% Load recorded response
Null_datafile = '/Volumes/Analysis/2015-03-09-2/data042-from-data038_streamed_nps';
neuronPath = [Null_datafile,sprintf('/data042-from-data038_streamed_nps.neurons')];
condDuration=10.6;
nConditions=6;
cond_str=[];



%% Do predictions - full

cell_list = datarun.cell_types{2}.cell_ids;

for cellID=cell_list
if(cellID==46)
    continue;
end
data1 = load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed/CellID_%d/Cell%d_full.mat',cellID,cellID));
data2 = load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed/CellID_%d/Cell%d_full_su4_5.mat',cellID,cellID));

extract_movie_response2
 if(~exist(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data042/OFF Parasol/null_pred_full_mix/cellID_%d',cellID),'dir'))
        mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data042/OFF Parasol/null_pred_full_mix/cellID_%d',cellID));
 end
print(hhh,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data042/OFF Parasol/null_pred_full_mix/cellID_%d/cellID_%d_mask.pdf',cellID,cellID));

fitGMLM_log=cell(5,1);
for isu=1:3
fitGMLM_log{isu} = data1.fitGMLM_full2_log{isu};
end
for isu=4:5
fitGMLM_log{isu} = data2.fitGMLM_full2_log{isu};
end
mask=data2.totalMaskAccept2;

for nSU=1:5
    fitGMLM=fitGMLM_log{nSU};
    fitGMLM = scramble_fitGMLM_mix(fitGMLM);
    for igain=[1,2,3,10]
[spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_2_light(cellID,nConditions,condDuration,cond_str,neuronPath);


% Make predictions
pred1=cell(nConditions,1);
for  icond=1:nConditions
movd = condMov{icond};
 maskedMov= filterMov(movd,mask,squeeze(ttf))*igain;
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
nTrials=30;
 pred1{icond}= predictGMLM_full(fitGMLM,maskedMov,nTrials)';
end


   % h=figure('Color','w','PaperSize',[10,7],'PaperPosition',[0 0 42 7]);
    close all
   hh=plot_record_prediction(spkCondColl,pred1);
% 
      
       print(hh,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data042/OFF Parasol/null_pred_full_mix/cellID_%d/cellID_%d_SU_%d_gain_%d.pdf',cellID,cellID,nSU,igain));
%   
save(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data042/OFF Parasol/null_pred_full_mix/cellID_%d/data_cellID_%d_SU_%d_gain_%d.mat',cellID,cellID,nSU,igain),'pred1','mask','ttf','spkCondColl','spkColl')
    end
end

end

%% Calculate PSTH and find correlation between recorded and predictions.


cell_list = datarun.cell_types{2}.cell_ids;

for cellID=cell_list

for nSU=1:5
    for igain=[1,2,3,10]
        clear spkCondColl pred1 
    load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data042/OFF Parasol/null_pred_full_mix/cellID_%d/data_cellID_%d_SU_%d_gain_%d.mat',cellID,cellID,nSU,igain));
    % predicted 
    convolve=150;
    len = 12720;
    binSz=1/1200;
    R2_log=[];
    PSTH_rec_log = cell(nConditions,1);
    PSTH_pred_log = cell(nConditions,1);
    
    for icond=1:nConditions
    [PSTH_pred,time]=calculate_psth_fcn2(convolve,binSz,len,pred1{icond});
    
    % Recorded 
    spkMat = makeSpikeMat(spkCondColl(icond).spksColl,binSz,len);
    [PSTH_rec,time]=calculate_psth_fcn2(convolve,binSz,len,spkMat);
     
   R2_log(icond) = R_2_value(PSTH_rec',PSTH_pred'); 
    PSTH_pred_log{icond} = PSTH_pred;
    PSTH_rec_log{icond}=PSTH_rec;
    end
    save(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data042/OFF Parasol/null_pred_full_mix/cellID_%d/correlation_cellID_%d_SU_%d_gain_%d.mat',cellID,cellID,nSU,igain),'R2_log','PSTH_pred_log','PSTH_rec_log')
    end
end
end

%% Make figures. - R2 value.

icond=1;
jcond=5;

icnt=0;
cellID_list= [];
for nSU=4%1:5
    for igain=3%[1,2,3,10]
        icnt=icnt+1;
        conda=[];condb=[];
        
        for cellID=cell_list
            aa = load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data042/OFF Parasol/null_pred_full/cellID_%d/correlation_cellID_%d_SU_%d_gain_%d.mat',cellID,cellID,nSU,igain));
            conda = [conda;aa.R2_log(icond)];
            condb=  [condb;aa.R2_log(jcond)];
            cellID_list = [cellID_list;cellID];
        end
        subplot(5,4,icnt);
        scatter(conda,condb,10,'filled');
        hold on;
        plot([0,1],[0,1],'g');
        axis square
        set(gca,'XTick',[0,1]);
        set(gca,'YTick',[0,1])
    end
end

%% Difference in original and null cell response difference triggered R-2 value

icond=1;
jcond=5;

icnt=0;
for nSU=1:5
    for igain=[1,2,3,10]
        icnt=icnt+1;
        conda=[];condb=[];
        
        for cellID=cell_list
            aa = load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data042/OFF Parasol/null_pred_full_scramble/cellID_%d/correlation_cellID_%d_SU_%d_gain_%d.mat',cellID,cellID,nSU,igain));
            
            p_rec1 = aa.PSTH_rec_log{icond};
            p_rec2 = aa.PSTH_rec_log{jcond};
            thr = prctile(abs(p_rec1-p_rec2),90);
            times =abs(p_rec1-p_rec2)>thr;
            
            
            conda = [conda;R_2_value(aa.PSTH_rec_log{icond}(times)',aa.PSTH_pred_log{icond}(times)')];
            condb=  [condb;R_2_value(aa.PSTH_rec_log{jcond}(times)',aa.PSTH_pred_log{jcond}(times)')];
        
        
        end
        subplot(5,4,icnt);
        scatter(conda,condb,10,'filled');
        hold on;
        plot([0,1],[0,1],'g');
        axis square
        set(gca,'XTick',[0,1]);
        set(gca,'YTick',[0,1])
    end
end
