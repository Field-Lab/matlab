  

WN_datafile = '2015-03-09-2/streamed/data038/data038';
Null_datafile = '/Volumes/Analysis/2015-03-09-2/data042-from-data038_streamed_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2015-03-09-2/data038/data038';
neuronPath = [Null_datafile,sprintf('/data042-from-data038_streamed_nps.neurons')];


datarun=load_data(WN_datafile)
datarun=load_params(datarun)

% 
 cellTypeId=[2]; % 1 for On Parasols, 2 for Off parasols
% InterestingCell_vis_id=[];
% for icellType=cellTypeId
%     icellTypecle
%     InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
% end

%ON parasol
%InterestingCell_vis_id = [421,857,962,1231,1592,1922,2611,3093,3121,3317,3590,3811,4187,4548,4727,5417,5568,5716,6091,6665,6737,6826,7473,7592];
%nullCond=3;

% OFF parasol
cellTypeId=2;
InterestingCell_vis_id = [842,1232,1233,1531,1981,2371,2596,2767,2806,3766,4052,4547,5326,5432,5551,5913,6287,6517,6694,6968,7172,7396,7518,7681]
nullCond = 5;

cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);
NullCells1=[2206,4548,6826,662]; % cell : 6826 ?
NullCells2=datarun.cell_types{1}.cell_ids;
NullCells3=[2596,1232,7172,842,4547,5911];
NullCells4=datarun.cell_types{2}.cell_ids;
NullCells5=[5768];

% InterestingCell_vis_id = [4548]%[662,6826,4548] %NullCells1;
condDuration=10.6;
nConditions=6;
GLM_fit_link= '/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data038/CellType_OFF parasol/';

LL_log=zeros(length(InterestingCell_vis_id),nConditions,2);
LL_log_glm=zeros(length(InterestingCell_vis_id),2,2);


for ref_cell_number=4%1:length(InterestingCell_vis_id); %11
    close all
    
    cellID=InterestingCell_vis_id(ref_cell_number)
    % Plot recorded raster
    [spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_2_light(cellID,nConditions,condDuration,cond_str,neuronPath);

    LL_log(ref_cell_number,:,:) = getLL_const_var_2015_03_09_2(spkCondColl,1);
    
     [resp1,resp2 ] = framework_null_expt_2015_03_09_2(cellID,GLM_fit_link,condMov{1}-0.5) % give movie between -0.5 to +0.5
     glm_resp(1).spksColl=resp1;glm_resp(2).spksColl = resp2;
     LL_log_glm(ref_cell_number,:,:) = getLL_const_var_2015_03_09_2(glm_resp,2);
     
     
end

LL_diff_orig = zeros(size(LL_log,1),2);
LL_diff_orig(:,1) = LL_log(:,1,2) - LL_log(:,1,1);
LL_diff_orig(:,2) = LL_log_glm(:,1,2) - LL_log_glm(:,1,1);



LL_diff_null= zeros(size(LL_log,1),2);
LL_diff_null(:,1) = LL_log(:,nullCond,2) - LL_log(:,nullCond,1);
LL_diff_null(:,2) = LL_log_glm(:,2,2) - LL_log_glm(:,2,1);

%save(sprintf('%sDiffLL.mat',GLM_fit_link),'LL_diff_orig','LL_diff_null','InterestingCell_vis_id','LL_log','LL_log_glm');

figure;plot(LL_log(:,1,2) - LL_log(:,1,1) , LL_log(:,3,2) - LL_log(:,3,1) ,'*');hold on;plot([0,0.8],[0,0.8],'g')

figure;plot(LL_diff_orig(:,2),LL_diff_orig(:,1),'*');hold on;plot(LL_diff_null(:,2),LL_diff_null(:,1),'*');hold on;plot([0,4],[0,4],'g');
xlabel('Framework');
ylabel('Data');
gmail('bhaishahster@gmail.com','Population Null update','Off Parasol done');

%% cell 1531 figure
close all
%[x1,y1] = plotSpikeRaster(resp2>0,'PlotType','vertline');
 [x1,y1] = plotSpikeRaster(spkCondColl(5).spksColl,'PlotType','vertline');
 x1=x1*120/20000;
h=figure('Color','w');
plot(x1/120,y1,'k');
xlim([0,1271/120]);
ylim([0,29]);
set(gca,'yTick',[]);
set(gca,'xTick',[]);
set(gca,'Visible','off')
set(h,'Position',pos)
s=hgexport('readstyle','singleRaster');
hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/Null_Figures_for_EJ/raster_framework_orig2.eps'),s);

%% 
onpar = load('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data038/CellType_ON parasol/DiffLL.mat');
offpar = load('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data038/CellType_OFF parasol/DiffLL.mat');


h= figure;scatter(onpar.LL_diff_null(:,2),onpar.LL_diff_null(:,1),30,1*ones(size(onpar.LL_diff_null,1),1)); hold on;

scatter(offpar.LL_diff_null(:,2),offpar.LL_diff_null(:,1),30,1*ones(size(onpar.LL_diff_null,1),1),'filled');hold on;
scatter(offpar.LL_diff_null(offpar.InterestingCell_vis_id==1531,2),offpar.LL_diff_null(offpar.InterestingCell_vis_id==1531,1),30,2,'filled');
%set(gca,'xTick',[]);
%set(gca,'yTick',[]);
hold on;plot([0,4],[0,4],'g');
xlabel('Framework')
ylabel('Data');
leg{1} ='On Parasol'; leg{2} ='Off Parasol';
legend(leg)
hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/Null_Figures_for_EJ/Population.eps'));

%%

% 
 cellTypeId=[1]; % 1 for On Parasols, 2 for Off parasols
% InterestingCell_vis_id=[];
% for icellType=cellTypeId
%     icellType
%     InterestingCell_vis_id=[InterestingCell_vis_id,datarun.cell_types{icellType}.cell_ids];
% end

%ON parasol
InterestingCell_vis_id = [421,857,962,1231,1592,1922,2611,3093,3121,3317,3590,3811,4187,4548,4727,5417,5568,5716,6091,6665,6737,6826,7473,7592];
nullCond=3;

% OFF parasol
% cellTypeId=2;
% InterestingCell_vis_id = [842,1232,1233,1531,1981,2371,2596,2767,2806,3766,4052,4547,5326,5432,5551,5913,6287,6517,6694,6968,7172,7396,7518,7681]
% nullCond = 5;

cellTypeUsed=cellTypeId*ones(length(InterestingCell_vis_id),1);
NullCells1=[2206,4548,6826,662]; % cell : 6826 ?
NullCells2=datarun.cell_types{1}.cell_ids;
NullCells3=[2596,1232,7172,842,4547,5911];
NullCells4=datarun.cell_types{2}.cell_ids;
NullCells5=[5768];

% InterestingCell_vis_id = [4548]%[662,6826,4548] %NullCells1;
condDuration=10.6;
nConditions=6;
GLM_fit_link= '/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data038/CellType_ON parasol/';

LL_log=zeros(length(InterestingCell_vis_id),nConditions,2);
LL_log_glm=zeros(length(InterestingCell_vis_id),2,2);


for ref_cell_number=1:length(InterestingCell_vis_id); %11
    close all
    
    cellID=InterestingCell_vis_id(ref_cell_number)
    % Plot recorded raster
    [spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_2_light(cellID,nConditions,condDuration,cond_str,neuronPath);

    LL_log(ref_cell_number,:,:) = getLL_const_var_2015_03_09_2(spkCondColl,1);
    
     [resp1,resp2 ] = framework_null_expt_2015_03_09_2(cellID,GLM_fit_link,condMov{1}-0.5) % give movie between -0.5 to +0.5
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
gmail('bhaishahster@gmail.com','Population Null update','On Parasol done');
