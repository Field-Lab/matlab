%% Kalmar non-linearity analysis 

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

%% Sub-unit NL

load('/Volumes/Lab/Users/bhaishahster/analyse_2015_10_29_2/d_add/PSTH_cond3_4.mat');
%figure;histogram(data(1).condba_rat(data(1).cells_select==1),'Normalization','probability');hold on;histogram(data(2).condba_rat(data(2).cells_select==1),'Normalization','probability');
dataRuns = dataRuns_OFF_additivity;
cellType=2;
ub = 0.5%0.5;
lb= 0.4%0.4;
% movies = movies_OFF_addivitiy;
cols = 'rb';

% load nls
if(cellType==2)
    load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_OFF.mat')
else
    load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/lnp_data002/data_nls2_ON.mat')
end


iidx = 1:length(data_nls);
mcellid = (data(cellType).condba_rat<=ub & data(cellType).condba_rat>=lb & data(cellType).cells_select'==1);
cids = data(cellType).cellIDs(mcellid);

nSU=4;
%figure('Color','w');
for icell=1:length(cids)
cids(icell)
Cell = load(sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/OFF Parasol_quad_MEL/Cell_%d.mat',cids(icell)));
fits  = Cell.fitGMLM_log{nSU};
for ifit=1:nSU
mn = mean(fits.data_act.kx{ifit});
std = sqrt(var(fits.data_act.kx{ifit}));
range = -2:0.01:2;
rr = range*std+mn;
plot(range,(rr.^2).*(rr>0),cols(cellType));
hold on;
end
pause(0.2)
end

%% output NL
% get response 

WN_datafile = '/Volumes/Analysis/2015-10-29-2/d00_36-norefit/data002/data002';


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);

cellID = 3348;%860;
Cell = load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/SUs_data002/OFF Parasol_quad_MEL/Cell_3348.mat');
nSU=4;
fits  = Cell.fitGMLM_log{nSU};
stim_length=1800;

lam = fits.data_act.lam;
% Load spikes
master_idx         = find(datarun.cell_ids == cellID);


% Spike loading
spikes=datarun.spikes{master_idx};

% make STA 3D ?
%glm_cellinfo.WN_STA = squeeze(sum(glm_cellinfo.WN_STA,3)); % Doubt!!!!!!!
clear cell_savename

% Align the spikes and the movies;
spikes_adj=spikes;
n_block=0;
for i=1:(length(datarun.triggers)-1)
    actual_t_start=datarun.triggers(i);
    supposed_t_start=n_block*100/120;
    idx1=spikes > actual_t_start;
    idx2=spikes < datarun.triggers(i+1);
    spikes_adj(find(idx2.*idx1))=spikes(find(idx2.*idx1))+supposed_t_start-actual_t_start;
    n_block=n_block+1;
end
clear spikes
spike.home=spikes_adj;
clear spikes_adj;


spksGen = zeros(stim_length*120,1);
for ispike=1:length(spike.home)
    spksGen(floor(spike.home(ispike)*120)+1)=1;
end
spksGen = spksGen(1:stim_length*120);

% plotNL 
dt=1/120;
th_list=[];
for iprc=0:0.5:100
 thr = prctile(lam,iprc);
th_list= [th_list;thr]; 
end

meanl=[];meanR=[];meanPredR = [];
for ith=1:length(th_list)-1
iidx = (lam>th_list(ith)) & (lam<=th_list(ith+1));
meanl = [meanl;mean(lam(iidx))/dt];
meanR = [meanR;mean(spksGen(iidx))/dt];
end

%%
% Cell2.meanl=meanl;Cell2.meanR=meanR;Cell2.th_list=th_list;

% Cell1.meanl=meanl;Cell1.meanR=meanR;Cell1.th_list=th_list;
%%

figure('Color','w');
loglog(Cell1.meanl,Cell1.meanR,'b');
hold on;
[X,N]=hist(Cell1.fitGMLM_log{nSU}.data_act.lam);
semilogx(N,10*X/sum(X))
%hold on;

loglog([Cell1.th_list(end-1),Cell1.th_list(end-1)]/dt,[0,35],'b');
hold on;
loglog(Cell2.meanl,Cell2.meanR,'r');hold on;
loglog([Cell2.th_list(end-1),Cell2.th_list(end-1)]/dt,[0,35],'r');

