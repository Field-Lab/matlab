addpath(genpath('../null_analyse/'));
addpath(genpath('../null_analyse/analyse_functions'));
%startup_null_analyse_tenessee
startup_null_analyse_bertha


%% Additivity experiment analysis 


nConditions=7;
condDuration=1200/4;
cond_str=cell(nConditions,1);
cond_str{1}='Original';
cond_str{2}='Null ';


dataRuns_OFF_additivity = [3,4,6,7,9,11,13];
dataRuns_ON_additivity = [3,5,6,8,10,12,13];
movies_OFF_addivitiy =[1,2,5,6,10,14,13];
movies_ON_additivity = [1,4,5,8,12,16,13];

location = '/Volumes/Analysis/2015-10-29-2/d00_36-norefit';

%% ASM fitting code

cellTypeId = 2%OFF:2 .. ON:1.. 
cellID_select= datarun.cell_types{cellTypeId}.cell_ids; % 51 ? [6741,2176,6961]

WN_datafile = '/Volumes/Analysis/2015-10-29-2/d00_36-norefit/data002/data002'
movie_xml = 'RGB-8-2-0.48-11111';
stim_length=1800;% 
cellID_list = [cellID_select];%1006,1411,2086,2341,2686,4096,4276,4831,5371,5566,5596,5866,7156,7216,7381,7490];
nSU_list= [1,2,3,4,5,7,10];
destination= 'pc2015_10_29_2_analysis_fits/SUs_data002/OFF '
sta_depth=30; % CAREFUL TO CHANGE IT TO LONG .. FOR SLOWER CELLS
contrast_factor=0.5;
[~] = get_gmlm_sta2(WN_datafile,movie_xml,stim_length,cellID_list,nSU_list,destination,sta_depth,contrast_factor)
save(['/Volumes/Lab/Users/bhaishahster/',destination,'/sus.mat'],'stas_t','stas_r');




%% ASM predictions with max lam code

%% compute op NL
WN_datafile = '/Volumes/Analysis/2015-10-29-2/d00_36-norefit/data002/data002';

datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun)
cellTypeId = 2;
InterestingCell_vis_id= datarun.cell_types{cellTypeId}.cell_ids; 

movie_xml = 'RGB-8-2-0.48-11111';
stim_length=1800;
cellIDs=InterestingCell_vis_id;
userSTA_depth=30;
destination = 'pc2015_10_29_2_analysis_fits/SUs_data002/OFF Parasol2';
contrast_factor=0.5;
save_location = 'pc2015_10_29_2_analysis_fits/SUs_data002/OFF Parasol_opnl_ARVO_nll2';

compute_op_NL_WN(WN_datafile,movie_xml,stim_length,cellIDs,userSTA_depth,destination,contrast_factor,save_location);

%%  ASM predictions with op NL code
   
    load('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_2/d_add/PSTH_cond3_4.mat');
    figure;histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');
    dataRuns = dataRuns_OFF_additivity;
    cellType=2;
    ub = 0.5;
    lb= 0.4;
    movies = movies_OFF_addivitiy;


% load nls
 if(cellType==2)
   load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_OFF.mat')
  else
    load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_ON.mat')
  end


  iidx = 1:length(data_nls);
  mcellid = iidx;%(data(cellType).condba_rat<=ub & data(cellType).condba_rat>=lb & data(cellType).cells_select'==1);
 cids = data(cellType).cellIDs(mcellid);
  
  
% load movies
interval=4;
condMov=cell(nConditions,1);
rawMovFrames=1200/(4);
icnt=0;
% make pixel histogram
for imov=movies
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-10-29-2/Visual/%d.rawMovie',imov),rawMovFrames,1);
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
cond_list=[3,4];
% iterate over different cells

sus_list = [1,2,3,4,5,6,7,8,9,10];
nSUs= length(sus_list);
pred_log = cell(length(mcellid),nSUs,length(cond_list));
R2_log = zeros(length(mcellid),nSUs,length(cond_list));
cellID_log=[];
nl_info_log=cell(length(mcellid),nSUs,length(cond_list));
error_cells=[];
metrics  = zeros(length(mcellid),nSUs,length(cond_list),4);

for imcell=1:length(mcellid)
    imcell
    
    display(sprintf('Starting: %d out of %d cells',imcell,length(mcellid)));
    current_cell=mcellid(imcell);
    cellID = cids(imcell);
    % load true responses .
    cellID_log=[cellID_log;cellID];
    if(data_nls(current_cell).cellID== cellID)
        display('cell IDs match');
    end
    cols='rkrkrkrkrkrkkrkrkrkr';
  %  try 
    condDuration=10;
    nConditions=1;
    for idata=1:length(dataRuns)
        Null_datafile = sprintf('%s/data0%02d',location,dataRuns(idata));
        neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
        [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    end
    
    close all
    
    if cellType==2
    cellData = load(sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/OFF Parasol2/Cell_%d',cellID));
     cellData_opnl =  load(sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/OFF Parasol_opnl_ARVO_nll/Cell_%d',cellID));
       
    else
    cellData = load(sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/ON Parasol2/Cell_%d',cellID));
  
    end
    
    for nSU=sus_list%1:nSUs
        nSU
        fitGMLM = cellData.fitGMLM_log{nSU};
        fitGMLM.NL = cellData_opnl.fit_nl(nSU).fit_nll;
        close all
        R2_pl=[];
        icnt=1;pred=cell(2,1);clear ss;nTrials=30;
        for icond=cond_list
            % real response
            realResp = makeSpikeMat(spkCondColl{icond}.spksColl,1/120,1200);
            convolve=7;
            binSz = 1/120;len=1200;
            [PSTH_rec,time]=calculate_psth_fcn2(convolve,binSz,len,realResp);
            ss(icnt)=spkCondColl{icond};
            [nSU,icond,sum(cellData.mask(:))]
            % predicted response
            [real_resp_bin1,time]=calculate_psth_fcn2(1,binSz,len,realResp);
            movd = (condMov{icond}-0.5);  % *2 if it is a low contrast movie!
            maskedMov= filterMov(movd,cellData.mask,squeeze(cellData.ttf));
            
            
            [pred{icnt},lam]= predictGMLM_bias_lr_opnl(fitGMLM,maskedMov,nTrials,1);pred{icnt} = pred{icnt}';
            %[pred{icnt},lam]= predictGMLM_bias_lr(fitGMLM,maskedMov,nTrials,1);pred{icnt} = pred{icnt}';
            %[pred{icnt},lam,nl_info_log{imcell,nSU,icnt}]= predictGMLM_bias_lr_optimal_op_nl(fitGMLM,maskedMov,nTrials,1,real_resp_bin1);
            %[pred{icnt},lam,nl_info_log{imcell,nSU,icnt}]= predictGMLM_gamma2_lr_optimal_op_nl(fitGMLM,maskedMov,nTrials,1,real_resp_bin1);
            %[pred{icnt},lam]= predictGMLM_gamma2_lr(fitGMLM,maskedMov,nTrials,2,1); pred{icnt}=pred{icnt}';
            %[pred{icnt},lam]= predictGMLM_full(fitGMLM,maskedMov,nTrials); pred{icnt}=pred{icnt}';
         
            nl_info_log{imcell,nSU,icnt} = 0;
            
            % [pred{icnt},lam,nl_info_log{imcell,nSU,icnt}]= predictGMLM_bias_lr(fitGMLM,maskedMov,nTrials,1);
            % pred{icnt} = pred{icnt}';
            % pause
           
            ss(icnt)=spkCondColl{icond};
            [PSTH_pred,time]=calculate_psth_fcn2(convolve,binSz,len,pred{icnt});
           
            R2_pl = [R2_pl;R_2_value(PSTH_rec(0.2*end:0.8*end)',PSTH_pred(0.2*end:0.8*end)')];
            PSTH_rec_log{icnt} = PSTH_rec;PSTH_pred_log{icnt}=PSTH_pred;
            
            R2_log(imcell,nSU,icnt)=R2_pl(end);
            
            metrics(imcell,nSU,icnt,:) = compute_raster_metrics(realResp,pred{icnt});
            p{1} =pred{icnt};
            pred_log(imcell,nSU,icnt) =p;
            icnt=icnt+1;
            
        end
        
        h=figure;%('Color','w','PaperSize',[10,7],'PaperPosition',[0 0 42 7]);
        PSTH_str.plot=1;PSTH_str.PSTH_rec_log = PSTH_rec_log;PSTH_str.PSTH_pred_log=PSTH_pred_log;
        %hh=plot_record_prediction(ss,pred,'PSTH_struct',PSTH_str);
        hh=plot_record_prediction(ss,pred);
        title(sprintf(' #SU:%d , %0.04f,%0.04f',nSU,R2_pl(1),R2_pl(2)));
        %title(sprintf(' %0.04f,%0.04f,%0.04f,%0.04f',R2_pl(1)));
        % xlim([8,10]);
         pause(1)
    end
%      catch
%          display(sprintf('Error in cell: %d, conditions: %d',imcell,icond));
%          error_cells = [error_cells,cids(imcell)];
%          
%      end
   save('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/OFF Parasol_opnl_ARVO/predictions_cond3_4_5_6_OFF_metrics_cutoff_axdybplx_negll.mat','metrics','pred_log','R2_log','error_cells','cellID_log','cond_list','sus_list');

 
end

%% see results
load('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_2/d_add/PSTH_cond3_4.mat');
figure;histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');
dataRuns = dataRuns_OFF_additivity;
cellType=2;%ub = 0.6; lb= 0.4;
lb=0.4;ub=0.6;
if(cellType==1)
    data_pred = load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002_optimalNL/predictions_cond3_4_5_6_ON_optimalNL_nSU10.mat')
else
    %data_pred = load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002_optimalNL/predictions_cond3_4_5_6_OFF_optimalNL_nSU10.mat')
    %data_pred = load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/predictions_cond3_4_5_6_OFF_quadraticSU.mat')
    %data_pred = load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002_optimalNL/predictions_cond3_4_5_6_OFF_metrics_optimalNL.mat')
     data_pred = load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/OFF Parasol_opnl_ARVO_nll/predictions_cond3_4_5_6_OFF_metrics_cutoff_axdybplx_negll.mat');
   
    %data_pred = load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/predictions_cond3_4_5_6_OFF_quadraticSU_metrics.mat')
end

iidx = 1:length(data(cellType).cellIDs);
cell_bin = (data(cellType).condba_rat<=ub & data(cellType).condba_rat>=lb & data(cellType).cells_select'==1);
mcellid = iidx(cell_bin);
cids = data(cellType).cellIDs(cell_bin);

figure;
icnt=0;
for isu=data_pred.sus_list
    icnt=icnt+1;
subplot(2,5,icnt)
plot(data_pred.R2_log(cell_bin,isu,1),data_pred.R2_log(cell_bin,isu,2),'.','markersize',5);
hold on;
plot([0,1],[0,1],'g')
xlim([0,1]);ylim([0,1]);
axis square
title(sprintf('nSU: %d ',isu));
end

figure;
for isu=data_pred.sus_list

    subplot(2,5,isu)
    for icond =[1,2]
    plot(data_pred.metrics(cell_bin,isu,icond,3),data_pred.metrics(cell_bin,isu,icond,1),'.');  
    hold on;
    end
    %legend('WN','Null','WN+Null var controlled','WN+Null var not controlled');
    title(sprintf('# SU %d:',isu));
end

% vary trials
clear vary_trials
icond=3;icnt=1;isu=5;

m_log=[]; 
for itrial=1:5
vary_trials(itrial).m_log=[];
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
    
    itrial=0;
    for nTrials = [4,6,8,12,20]
        itrial=itrial+1;
        realResp = makeSpikeMat(spkCondColl{icond}.spksColl,1/120,1200);
        pred=data_pred.pred_log(icell,isu,icnt);
        if(~isempty(pred{1}))
        m = compute_raster_metrics(realResp(1:floor(30/nTrials):nTrials*floor(30/nTrials),:),pred{1}(1:floor(30/nTrials):nTrials*floor(30/nTrials),:));
        %m_log=[m_log;[m(3),m(1)]];
        else
            m=[0,0,0,0];
        end
        vary_trials(itrial).m_log=[vary_trials(itrial).m_log;[m(3),m(1)]];
    end
    
end

figure;
for itrial = 1:length(vary_trials)
    hold on;
    plot(vary_trials(itrial).m_log(:,1),vary_trials(itrial).m_log(:,2),'k.','markersize',15);
end

for icnt =[1,2]
    plot(data_pred.metrics(cell_bin,isu,icnt,3),data_pred.metrics(cell_bin,isu,icnt,1),'.','markersize',15,'Color',cols(icnt,:));
    hold on;
end


%legend('WN','Null','WN+Null var controlled','WN+Null var not controlled');
title(sprintf('# SU %d:',isu));


% be more selective about which cells to plot!
% cell_bin2 = cell_bin'  &  sum(data_pred.metrics(:,:,1,1)~=0,2)>0 & sum(data_pred.metrics(:,:,1,3)~=0,2)>0 &  sum(data_pred.metrics(:,:,2,1)~=0,2)>0 & sum(data_pred.metrics(:,:,2,3)~=0,2)>0 ;  

% convert cellids_select to cell_bin2
load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/OFF Parasol_select.mat')
cell_bin2=zeros(numel(data(cellType).cellIDs),1);
for icell=1:numel(data(cellType).cellIDs)
if(sum(data(cellType).cellIDs(icell) == cids_select)>0)
cell_bin2(icell)=1;
end
end
cell_bin2 = logical(cell_bin2);

figure;
   % cols = distinguishable_colors(6);
cols=[0,0,0;1,0,0];
   for icond=1:2;

    for iisu=1:length(data_pred.sus_list(1:end-1))
        isu = data_pred.sus_list(iisu);
        jsu = data_pred.sus_list(iisu+1);
        metric_ori = data_pred.metrics(cell_bin2,isu,icond,1)./data_pred.metrics(cell_bin2,isu,icond,3);
        metric_end = data_pred.metrics(cell_bin2,jsu,icond,1)./data_pred.metrics(cell_bin2,jsu,icond,3);
        for icell=1:sum(cell_bin2)
            a=plot([isu+(icond-1)*0.1,jsu+(icond-1)*0.1],[metric_ori(icell),metric_end(icell)],'-*','Color',cols(icond,:));
            a.Color(4) = 0.25;
            hold on;
            
%             if(iisu==4)
%             text(4,metric_ori(icell),num2str(cids_select(icell)));
%             end
        end
        
    end
    
    
end

for icond=1:2;
    for iisu=1:length(data_pred.sus_list(1:end-1))
        isu = data_pred.sus_list(iisu);
        jsu = data_pred.sus_list(iisu+1);
        metric_ori = data_pred.metrics(cell_bin2,isu,icond,1)./data_pred.metrics(cell_bin2,isu,icond,3);
        metric_end = data_pred.metrics(cell_bin2,jsu,icond,1)./data_pred.metrics(cell_bin2,jsu,icond,3);
        plot([isu+(icond-1)*0.1,jsu+(icond-1)*0.1],[mean(metric_ori(~isnan(metric_ori))),mean(metric_end(~isnan(metric_end)))],'-*','Color',cols(icond,:),'LineWidth',3);
        hold on;
    end
end
xlim([0.5,10.5]);

%% plot rasters


   for  cellID =6275%cids_select(randperm(length(cids_select)))';
       cellID
condDuration=10;
    nConditions=1;
    cond_str=[];
    
    for idata=1:length(dataRuns)
    Null_datafile = sprintf('/Volumes/Analysis/2015-10-29-2/d00_36-norefit/data0%02d',dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
    [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    close all
    end
    
    % WN
    susplot=[1,3,5];predCell = cell(length(susplot),1);
    clear ss
    for nsuplot=1:length(susplot)
    ss(nsuplot)=spkCondColl{3};
    predCell{nsuplot} = data_pred.pred_log{data_pred.cellID_log==cellID,susplot(nsuplot),1};
    end
    figure;
    hh=plot_record_prediction(ss,predCell);
    title(sprintf('%f,%f,%f',data_pred.metrics(data_pred.cellID_log==cellID,susplot(1),1,1)./data_pred.metrics(data_pred.cellID_log==cellID,susplot(1),1,3),data_pred.metrics(data_pred.cellID_log==cellID,susplot(2),1,1)./data_pred.metrics(data_pred.cellID_log==cellID,susplot(2),1,3),data_pred.metrics(data_pred.cellID_log==cellID,susplot(3),1,1)./data_pred.metrics(data_pred.cellID_log==cellID,susplot(3),1,3)));
      xlim([3,9])
      
    % Null
    predCell = cell(length(susplot),1);
    clear ss
    for nsuplot=1:length(susplot)
    ss(nsuplot)=spkCondColl{4};
    predCell{nsuplot} = data_pred.pred_log{data_pred.cellID_log==cellID,susplot(nsuplot),2};
    end
    figure;
    hh=plot_record_prediction(ss,predCell);
    title(sprintf('%f,%f,%f',data_pred.metrics(data_pred.cellID_log==cellID,susplot(1),2,1)./data_pred.metrics(data_pred.cellID_log==cellID,susplot(1),2,3),data_pred.metrics(data_pred.cellID_log==cellID,susplot(2),2,1)./data_pred.metrics(data_pred.cellID_log==cellID,susplot(2),2,3),data_pred.metrics(data_pred.cellID_log==cellID,susplot(3),2,1)./data_pred.metrics(data_pred.cellID_log==cellID,susplot(3),2,3)));
   xlim([3,9])
   pause
   end