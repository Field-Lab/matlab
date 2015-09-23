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

 cellID=1531

data1 = load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed/CellID_%d/Cell%d_full.mat',cellID,cellID));
data2 = load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed/CellID_%d/Cell%d_full_su4_5.mat',cellID,cellID));
user_STA_depth=30;
extract_movie_response2
%  if(~exist(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data042/OFF Parasol/null_pred_full_mix/cellID_%d',cellID),'dir'))
%         mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data042/OFF Parasol/null_pred_full_mix/cellID_%d',cellID));
%  end
% print(hhh,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data042/OFF Parasol/null_pred_full_mix/cellID_%d/cellID_%d_mask.pdf',cellID,cellID));

fitGMLM_log=cell(5,1);
for isu=1:3
fitGMLM_log{isu} = data1.fitGMLM_full2_log{isu};
end
for isu=4:5
fitGMLM_log{isu} = data2.fitGMLM_full2_log{isu};
end
mask=data2.totalMaskAccept2;
%%
 nSU=4
    fitGMLM=fitGMLM_log{nSU};
   % fitGMLM = scramble_fitGMLM_mix(fitGMLM);

   
[spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_2_light(cellID,nConditions,condDuration,cond_str,neuronPath);
   fr=[];
   gain_range = [1:0.2:5];
    for igain=gain_range

close all;
% Make predictions
pred1=cell(1,1);
icond=5%1:nConditions
movd = condMov{icond};
 maskedMov= filterMov(movd,mask,squeeze(ttf))*igain;
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
nTrials=29;
 [pred1{1},lam]= predictGMLM_full(fitGMLM,maskedMov,nTrials);
 pred1{1} = pred1{1}';
 fr =[fr;sum(pred1{1}(:))/(size(pred1{1},1) * size(pred1{1},2)/1200)];
 
    end
figure;
plot(gain_range,fr);
hold on;
plot(gain_range,spkCondColl(icond).avgSpkRate *ones(size(gain_range)));

realResp = makeSpikeMat(spkCondColl(icond).spksColl,1/120,1272);
fr= mean(realResp,1)/(1/120);

igain=input('What gain to choose?');
pred1=cell(1,1);
icond=5%1:nConditions
movd = condMov{icond};
 maskedMov= filterMov(movd,mask,squeeze(ttf))*igain;
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
nTrials=29;
 [pred1{1},lam]= predictGMLM_full(fitGMLM,maskedMov,nTrials);
 pred1{1} = pred1{1}';
 
   % h=figure('Color','w','PaperSize',[10,7],'PaperPosition',[0 0 42 7]);
    close all
    ss = spkCondColl(icond);
   hh=plot_record_prediction(ss,pred1);
% 
      
 %      print(hh,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data042/OFF Parasol/null_pred_full_mix/cellID_%d/cellID_%d_SU_%d_gain_%d.pdf',cellID,cellID,nSU,igain));
%   
%save(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data042/OFF Parasol/null_pred_full_mix/cellID_%d/data_cellID_%d_SU_%d_gain_%d.mat',cellID,cellID,nSU,igain),'pred1','mask','ttf','spkCondColl','spkColl')
    

