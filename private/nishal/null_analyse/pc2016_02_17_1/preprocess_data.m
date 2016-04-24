%% pc2016-02-17-1

% Condition strings
nConditions=3;
condDuration=2400/4;
cond_str=cell(nConditions,1);
cond_str{1}='Original';
cond_str{2}='Null ';

interestingConditions=[1,2,3];

dataRuns_OFF_additivity = [31,32,34,38];
dataRuns_ON_additivity = [31,33,35,38];


WN_datafile = '/Volumes/Analysis/2016-02-17-1/d28_51-norefit/data028/data028';
WN_datafile_streamed = '/Volumes/Analysis/2016-02-17-1/streamed/data028/data028';
location = '/Volumes/Analysis/2016-02-17-1/d28_51-norefit';

%%


dataRuns =dataRuns_OFF_additivity;
cellTypeId = 2;

datarun_s=load_data(WN_datafile_streamed);
datarun_s=load_params(datarun_s);
datarun_s = load_sta(datarun_s);
InterestingCell_vis_id_streamed=datarun_s.cell_types{cellTypeId}.cell_ids; 


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun = load_sta(datarun);

InterestingCell_vis_id=datarun.cell_types{cellTypeId}.cell_ids; 
cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);

condDuration=20;
nConditions=1;
conda=1;condb=2;

cols='rkrkrkrkrkrkkrkrkrkr';
spkCondColl=cell(length(dataRuns),1);
conda_var=[];condb_var=[];conda_mean=[];condb_mean=[];cells_select=[];
for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    cellID= InterestingCell_vis_id(ref_cell_number)
        
    for idata=1:length(dataRuns)
    Null_datafile = sprintf('%s/data0%02d',location,dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
   % [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_12_18_2_light_photons(cellID,nConditions,condDuration,cond_str,neuronPath);
    
    end
    
    
    h=figure;
    for irun = 1:length(dataRuns)
    plot(spkCondColl{irun}.xPoints/20000,spkCondColl{irun}.yPoints - (irun-1)*30,cols(irun));
    hold on;
    end
    set(gca,'yTick',[]);
    ylim([-(length(dataRuns)-1)*30,30]);
    InterestingCell_vis_id(ref_cell_number)
    
    % compare structure in rasters in conditions a and b
      convolve=150;
    binSz = 1/1200;len=12000;
    realResp = makeSpikeMat(spkCondColl{conda}.spksColl, binSz,len);
    [PSTH_reca,time]=calculate_psth_fcn2(convolve,binSz,len,realResp);
        
    realResp = makeSpikeMat(spkCondColl{condb}.spksColl,binSz,len);
    [PSTH_recb,time]=calculate_psth_fcn2(convolve,binSz,len,realResp);
    
    
    conda_var = [conda_var,sqrt(var(PSTH_reca(0.25*end:end)))];
    condb_var = [condb_var,sqrt(var(PSTH_recb(0.25*end:end)))];
    
    conda_mean = [conda_mean,sqrt(mean(PSTH_reca(0.25*end:end)))];
    condb_mean = [condb_mean,sqrt(mean(PSTH_reca(0.25*end:end)))];
     
    close all
    
    
   
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(6,3,[1:9]);
    for irun = 1:length(dataRuns)
    plot(spkCondColl{irun}.xPoints/20000,spkCondColl{irun}.yPoints - (irun-1)*30,cols(irun));
    hold on;
    end
    set(gca,'yTick',[]);
    ylim([-(length(dataRuns)-1)*30,30]);
    
    subplot(6,3,[10,11,12]);
   plot(time,PSTH_reca);hold on;plot(time,PSTH_recb);
    title(sprintf('%f,%f',conda_var(end),condb_var(end)));

    subplot(6,3,[13,16]);
    imagesc(mean(datarun.stas.stas{datarun.cell_ids==cellID}(:,:,:,27),3));
    colormap gray
    axis image;set(gca,'xTick',[]);set(gca,'yTick',[]);
   
    subplot(6,3,[14,17]);
   plot_rf_fit_nishal(datarun, InterestingCell_vis_id,'magnify',20,'fill_color',[0,1,1],'alpha',0.3,'fill',false,'edge',true,'labels',false,'edge_color',[1,0,0]);
   hold on;
    plot_rf_fit_nishal(datarun, cellID,'magnify',20,'fill_color',[0,1,1],'alpha',0.3,'fill',true,'edge',true,'labels',false,'edge_color',[1,0,0]);
   axis image;set(gca,'xTick',[]);set(gca,'yTick',[]);
  title('Concat class.');
      
    subplot(6,3,[15,18]);
   plot_rf_fit_nishal(datarun_s, InterestingCell_vis_id_streamed,'magnify',20,'fill_color',[0,1,1],'alpha',0.3,'fill',false,'edge',true,'labels',false,'edge_color',[1,0,0]);
   axis image;set(gca,'xTick',[]);set(gca,'yTick',[]);
 title('Streamed class.');
 
   suptitle(sprintf('CellID: %d',cellID));
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


cells_select = [cells_select;cell_ch];
end

condba_rat = condb_var./conda_var;

data(cellTypeId).conda_var=conda_var;
data(cellTypeId).condb_var=condb_var;
data(cellTypeId).condba_rat=condba_rat;
data(cellTypeId).cells_select = cells_select;
data(cellTypeId).conda_mean =conda_mean;
data(cellTypeId).condb_mean =condb_mean;
data(cellTypeId).cellIDs = InterestingCell_vis_id;
%%
figure('Color','w');
 histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability','NumBins',20);hold on;
 histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability','NumBins',20);
%xlim([0,1]);
legend('OFF parasol','ON parasol');
xlabel('strcut in Null/ struct in WN');

%% 
figure('Color','w');
scatter(data(2).conda_var(data(2).cells_select==1) ,data(2).condb_var(data(2).cells_select==1),'filled');
hold on;
scatter(data(1).conda_var(data(1).cells_select==1) ,data(1).condb_var(data(1).cells_select==1),'filled');
hold on;
plot([0,0.018],[0,0.018],'g');
axis equal
legend('OFF parasol','ON parasol','equality');
xlabel('structure in WN');
ylabel('structure in null');