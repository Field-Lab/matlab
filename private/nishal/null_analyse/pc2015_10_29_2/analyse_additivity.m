%% Initial stuff
%%
% Condition strings
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
location = '/Volumes/Analysis/2015-10-29-2/d00_36-norefit';



WN_datafile = '/Volumes/Analysis/2015-10-29-2/d00_36-norefit/data001/data001';


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);
    %% prediction using sub-unit model
cellID = 3348;
movies=movies_OFF_addivitiy;  
dataRuns = dataRuns_OFF_additivity;

% load movies
interval=4;
condMov=cell(nConditions,1);
rawMovFrames=1200/(4);
icnt=0;
% make pixel histogram
for imov=movies
    [stim,height,width,header_size] = get_raw_movie([path,sprintf('/stimulus/%d.rawMovie',imov)],rawMovFrames,1);
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

% load true responses .

location = [path,'/analysis_data'];%'/Volumes/Analysis/2015-10-29-2/d00_36-norefit';

cols='rkrkrkrkrkrkkrkrkrkr';
       
condDuration=10;
nConditions=1;
    for idata=1:length(dataRuns)
    Null_datafile = sprintf('%s/data0%02d',location,dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
    [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    end
    
location = [path,'/fits/'];%'/Volumes/Analysis/2015-10-29-2/d00_36-norefit';
       cellData = load(sprintf('%s/Cell_%d',location,cellID));
   for    nSU=4%1:10;
    fitGMLM = cellData.fitGMLM_log{nSU};
 
    

    R2_pl=[];
    icnt=1;pred=cell(2,1);clear ss;nTrials=30;
    for icond=[1,3,4]
        
         movd = (condMov{icond}-0.5);  % *2 if it is a low contrast movie!
         maskedMov= filterMov(movd,cellData.mask,squeeze(cellData.ttf));
         %  maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
        [pred{icnt},lam,kx]= predictGMLM_bias_lr(fitGMLM,maskedMov,nTrials,1);
        pred{icnt}=pred{icnt}';
%        pred{icnt} = %response_lnp(data_nls(current_cell),gain*(condMov{icond}-0.5),nTrials);    
        ss(icnt)=spkCondColl{icond};
        
        realResp = makeSpikeMat(spkCondColl{icond}.spksColl,1/120,1200);
        convolve=15;
        binSz = 1/120;len=1200;
        [PSTH_rec,time]=calculate_psth_fcn2(convolve,binSz,len,realResp);
        pred{icnt}(pred{icnt}>100000)=100000;
        [PSTH_pred,time]=calculate_psth_fcn2(convolve,binSz,len,pred{icnt});

        R2_pl = [R2_pl;R_2_value(PSTH_rec(0.2*end:0.8*end)',PSTH_pred(0.2*end:0.8*end)')];
        
        %% plot sub-unit NL
        
        th=0.3;
        figure('Color','w');
        subplot(2,1,1);
        plot(PSTH_rec);hold on;plot(PSTH_rec*0+ th*max(PSTH_rec),'g');
        
        subplot(2,1,2);
        for ifit=1:nSU
            mn = mean(kx{ifit});
            std = sqrt(var(kx{ifit}));
            range = -3:0.01:3;
            %plot(range,exp(range*std+mn),'b');
            hold on;
            
            %when high firing events, whats the generator signal? 
            signal= kx{ifit}(PSTH_rec>th*max(PSTH_rec));
            [x,n] = hist((signal-mn)/std);
            hold on;%plot(n,x/sum(x),'r');
            plotyy(n,x/sum(x),range,exp(range*std+mn));
            
        end
        
        hold on;
        plotyy(range,normpdf(range,0,1),0,0);
        
        
        
        %%
         icnt=icnt+1;
        
    end
    
h=figure('Color','w','PaperSize',[10,7],'PaperPosition',[0 0 42 7]);
PSTH_str.plot=1;PSTH_str.PSTH_rec_log = PSTH_rec_log;PSTH_str.PSTH_pred_log=PSTH_pred_log;
hh=plot_record_prediction(ss,pred);
title(sprintf('nSU : %d',nSU));

print(h,sprintf('/Users/bhaishahster/Google Drive/new york/nSU_%d.pdf',nSU));
%title(sprintf('gain: %d, %0.04f,%0.04f',gain,R2_pl(1),R2_pl(2)));
% xlim([8,10]);
% pause


   end
  

%save('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/predictions_cond3_4_OFF.mat','pred_log','R2_log','gain_list','cellID_log');
 
%% compare WN and WN + Null stimuli cell responses


load('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_2/d_add/PSTH_cond3_4.mat');
figure;histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');
dataRuns = dataRuns_ON_additivity;
cellType=1;%ub = 0.6; lb= 0.4;
%lb=0.4;ub=0.6; % OFF
lb=0.15;ub=0.3; % ON

iidx = 1:length(data(cellType).cellIDs);
cell_bin = (data(cellType).condba_rat<=ub & data(cellType).condba_rat>=lb & data(cellType).cells_select'==1);
mcellid = iidx(cell_bin);
cids = data(cellType).cellIDs(cell_bin);

% vary conditions
clear var_cond
icond=3;
jcond_list = [4,5,6,7];
for jcond = jcond_list
var_cond(jcond).m_log=[];
end

for icell = mcellid
    cellID = data(cellType).cellIDs(icell);
    cellID
    condDuration=10;
    nConditions=1;
    close all
    for idata=1:length(dataRuns)
        Null_datafile = sprintf('%s/data0%02d',location,dataRuns(idata));
        neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
        [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    end    
    

for jcond=jcond_list
    
    realResp1 = makeSpikeMat(spkCondColl{icond}.spksColl,1/120,1200);
    realResp2 = makeSpikeMat(spkCondColl{jcond}.spksColl,1/120,1200);
    
    m = compute_raster_metrics(realResp1,realResp2);
    var_cond(jcond).m_log=[var_cond(jcond).m_log;[m(1),m(2),m(3),m(4)]];
end

end

figure;
for jcond=jcond_list
   plot(var_cond(jcond).m_log(:,3),var_cond(jcond).m_log(:,1),'.');
   hold on;
end
legend('Null','WN+Null var controlled','WN+Null var not controlled','WN again')

figure;
for jcond=jcond_list
    hold on;
  histogram(var_cond(jcond).m_log(:,1));
   hold on;
end
legend('Null','WN+Null var controlled','WN+Null var not controlled','WN again')

