%% Kalmar non-linearity analysis - new york

%% Sub-unit NL
path = '~/Google Drive/new york'
Cell1 = load([path,'/fits/Cell_3348.mat']);
Cell2 = load([path,'/fits/Cell_860.mat']);

nSU=4;
fits  = Cell1.fitGMLM_log{nSU};
figure('Color','w');
for ifit=1:nSU
mn = mean(fits.data_act.kx{ifit});
std = sqrt(var(fits.data_act.kx{ifit}));
range = -2:0.01:2;
plot(range,(range*std+mn),'b');
hold on;
end

fits  = Cell2.fitGMLM_log{nSU};
for ifit=1:nSU
mn = mean(fits.data_act.kx{ifit});
std = sqrt(var(fits.data_act.kx{ifit}));
range = -2:0.01:2;
plot(range,(range*std+mn),'r');
hold on;
end
hold on;
plot(range,10*normpdf(range,0,1),'k');

%% output NL
% get response 

javaaddpath('/Users/bhaishahster/Dropbox/Lab/Applications/Vision.app/Contents/Resources/Java/Vision.jar');
WN_datafile = ['/Users/bhaishahster/Google Drive/new york/analysis_data/data002/data002'];


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);

cellID = 3348;%860;
Cell = load([path,'/fits/Cell_3348.mat']);
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
plot(Cell1.meanl,Cell1.meanR,'b');
hold on;
plot([Cell1.th_list(end-1),Cell1.th_list(end-1)]/dt,[0,35],'b');
hold on;
plot(Cell2.meanl,Cell2.meanR,'r');hold on;
plot([Cell2.th_list(end-1),Cell2.th_list(end-1)]/dt,[0,35],'r');
