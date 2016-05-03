    %% prediction using optimal output non-linearity to match the response
    
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
        %data_pred = load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/predictions_cond3_4_5_6_OFF_quadraticSU_metrics.mat')
        data_pred = load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/OFF Parasol_opnl_ARVO_nll/predictions_cond3_4_5_6_OFF_metrics_cutoff_axdybplx_negll.mat');
   
    end

iidx = 1:length(data(cellType).cellIDs);
cell_bin = (data(cellType).condba_rat<=ub & data(cellType).condba_rat>=lb & data(cellType).cells_select'==1);
cids = data(cellType).cellIDs(cell_bin);
  


% load nls
%  if(cellType==2)
%    load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_OFF.mat')
%  else
%     load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_ON.mat')
%  end
% 

 
    % load sta 
    datafile = '/Volumes/Analysis/2015-10-29-2/d00_36-norefit/data002/data002.params';
    datarun=load_data(datafile)
    datarun=load_sta(datarun)
    datarun=load_params(datarun)
    
    datafile2 = '/Volumes/Analysis/2015-10-29-2/streamed/data001/data001.params';
    datarun2=load_data(datafile2)
    datarun2=load_sta(datarun2)
    datarun2=load_params(datarun2)
    %%
    cids_select=[];
for cellID=4475%cids
    %pause
       % plot rasters
    condDuration=10;
    nConditions=1;
    cond_str=[];
    
    for idata=1:length(dataRuns)
    Null_datafile = sprintf('/Volumes/Analysis/2015-10-29-2/d00_36-norefit/data0%02d',dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
    [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    close(h)
    end
    
    close all
    figure('units','normalized','outerposition',[0 0 1 1]);
    mask= load(sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/OFF Parasol2/Cell_%d',cellID),'mask');
    mask = mask.mask;
    [r,c] = find(mask==1);
    
    %get_cell_ids(datarun,type_name) % Vision IDs - used for loading fitted
    %STAs!
    stas=datarun.stas.stas(datarun.cell_ids==cellID);
    stas=stas{1};
    
    xx = double(squeeze(mean(stas(:,:,:,28),3)));
    subplot(5,3,1);
    imagesc(xx);colormap gray;axis image;
    for ir =1:length(r)
        hold on;
        plot(r(ir),c(ir),'r.','MarkerSize',5);
    end
    set(gca,'visible','off');
    caxis([min(stas(:)),max(stas(:))]);
    xlim([min(r)-2,max(r)+2]);
    ylim([min(c)-2,max(c)+2]);
    %title('data002 and significant stixels!');

    
     subplot(5,3,2);
    imagesc(xx);colormap gray;axis image;
    for ir =1:length(r)
        hold on;
        plot(r(ir),c(ir),'r.','MarkerSize',1);
    end
    caxis([min(stas(:)),max(stas(:))]);
    set(gca,'visible','off');
    suptitle(sprintf('Cell ID: %d',cellID));

%     
%     cellID_streamed =input(sprintf('cell ID in streamed/data001 for cell %d in concatenated analysis?',cellID));
%     stas=datarun2.stas.stas(datarun2.cell_ids==cellID_streamed);icell=1;
%     st_temp=zeros(size(stas{1},2),size(stas{1},1),1,size(stas{1},4)); % DOUBT .. Could be a reason for things to fail!!!!!
%     for itime=1:size(stas{1},4)
%         st_temp(:,:,:,itime)=mean(stas{1}(:,:,:,end-itime+1),3)';
%     end
%     stas_new{icell}=st_temp;
% 
%     stass=stas{1};
%     cell_params.STAlen=14;
%     [stas_clipped,totalMaskAccept,CellMasks]= clipSTAs(stas_new,cell_params);
%     [r,c] = find(totalMaskAccept==1);
%     xx = double(squeeze(mean(stass(:,:,:,26),3)));
%     subplot(1,2,2);
%     imagesc(xx);colormap gray;axis image;
%     for ir =1:length(r)
%         hold on;
%         plot(r(ir),c(ir),'r.','MarkerSize',20);
%     end
%     xlim([min(r)-2,max(r)+2]);
%     ylim([min(c)-2,max(c)+2]);
%     
 
    % WN
    susplot=[2,5,7];predCell = cell(length(susplot),1);
    clear ss
    for nsuplot=1:length(susplot)
    ss(nsuplot)=spkCondColl{3};
    predCell{nsuplot} = data_pred.pred_log{data_pred.cellID_log==cellID,susplot(nsuplot),1};
    end
    
    subplot(5,3,[4,5,6,7,8,9]);
    hh=plot_record_prediction(ss,predCell);
    title(sprintf('%f,%f,%f',data_pred.R2_log(data_pred.cellID_log==cellID,susplot(1),1),data_pred.R2_log(data_pred.cellID_log==cellID,susplot(2),1),data_pred.R2_log(data_pred.cellID_log==cellID,susplot(3),1)));
    
    % Null
    susplot=[1,4,7];predCell = cell(length(susplot),1);
    clear ss
    for nsuplot=1:length(susplot)
    ss(nsuplot)=spkCondColl{4};
    predCell{nsuplot} = data_pred.pred_log{data_pred.cellID_log==cellID,susplot(nsuplot),2};
    end
    subplot(5,3,[10,11,12,13,14,15]);
    hh=plot_record_prediction(ss,predCell);
    title(sprintf('%f,%f,%f',data_pred.R2_log(data_pred.cellID_log==cellID,susplot(1),2),data_pred.R2_log(data_pred.cellID_log==cellID,susplot(2),2),data_pred.R2_log(data_pred.cellID_log==cellID,susplot(3),2)));
   
   
    
    % dialogue if a good cell or a bad cell
% Construct a questdlg with two options
cell_ch=1;
choice = questdlg('Is it a good cell?', ...
	'Cell choice', ...
	'good cell','bad cell','good cell');
% Handle response
switch choice
    case 'good cell'
        disp([choice ' selected'])
        cell_ch = 1;
    case 'bad cell'
        disp([choice ' selected'])
        cell_ch = 0;
end
if(cell_ch==1)
cids_select=[cids_select;cellID];
end

end
