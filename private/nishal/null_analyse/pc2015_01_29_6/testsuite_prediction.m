

%% Load STA
opt=struct('load_all',true);
datarun=load_data(WN_datafile,opt)
datarun=load_sta(datarun)
datarun=load_params(datarun)
datarun=load_ei(datarun,'all','array_type',519);
%get_cell_ids(datarun,type_name) % Vision IDs - used for loading fitted
%STAs!

vision_id=InterestingCell_vis_id(ref_cell_number);
idx=[1:length(datarun.cell_ids)];
matlab_id=idx(datarun.cell_ids==vision_id);
cell_ei=datarun.ei.eis{matlab_id};
ssta=datarun.stas.stas{matlab_id};

%% Fit STA
addpath(genpath('~/GITs/matlab/code')); 
addpath(genpath('~/GITs/matlab/private/nishal/fwdfittingfunctions'));
fit_info = fit_sta(ssta)
%datarun = compute_sta_fits(datarun, 'all');

%params=datarun.matlab.sta_fits{matlab_id}.initial_params;
params=fit_info.initial_params;
full_fit = sta_fit_function(params);

figure;
for itime=27%30:-1:1
    itime
    subplot(2,1,1);
    imagesc(squeeze(sum(ssta(:,:,:,itime),3)));
    colormap gray
    axis image
    colorbar
    caxis([min(ssta(:)), max(ssta(:))]);

    
    subplot(2,1,2);
    imagesc(squeeze(sum(full_fit(:,:,:,itime),3)));
    colormap gray
    axis image
    colorbar 
    caxis([min(ssta(:)), max(ssta(:))]);
 
 %    pause;
    
end

%%

% NOTE: cell_types{i}.name to be chosen to include the target cell of interest. 
addpath(genpath('../../../../code/'));
datarun = get_sta_summaries(datarun, datarun.cell_types{cellTypeUsed(ref_cell_number)}.name, ...
    'verbose',0,'keep_stas',0,'keep_rfs',1,'fig_or_axes',[],...
    'marks_params',struct( ...
        'strength','vector length', 'filter', fspecial('gauss',15,0.7), ...
        'thresh',5,'robust_std_method',1));
movie_spec=WN_mov; %RGB-8-1-0.48-11111.xml';%RGB-8-1-0.48-11111-80x40
datarun = load_java_movie(datarun, movie_spec);
%datarun = get_snls(datarun, datarun.cell_ids(get_cell_indices(datarun, datarun.cell_types{1}.name)),'frames',-18:0,'stimuli',[],'new',true);
datarun = get_snls(datarun,vision_id,'frames',-18:0,'stimuli',[],'new',true);
% NOTE: TODO : have to Vary significant stixels, STA length (in frames parameter which is set at -18 right now !)and see how shape of non-linearity changes!

cellID=matlab_id;
gen=datarun.stas.snls{cellID}.gen_signal;
spks=datarun.stas.snls{cellID}.spikes;

figure;
subplot(2,2,1);
hist(spks);
subplot(2,2,2);
scatter(gen,spks);
hold on;
x=[-1:0.01:1];
N=@(x) exp(datarun.stas.snls{cellID}.fit_params.a*x +datarun.stas.snls{cellID}.fit_params.b);
plot(x,N(x),'r');

subplot(2,2,4);
hist(gen,100);

% int=[-1:0.01:1];
% use gen and make non-uniform grid
int = quantile(gen,100);

gen_log=[];
p_spk=[];
error_bar=[];
val_ids=(gen<=int(1));
gen_log=[gen_log;mean(gen(val_ids))];
p_spk=[p_spk;mean(spks(val_ids))];
error_bar=[error_bar;sqrt(var(spks(val_ids)))/sqrt(sum(val_ids))];

for idx=1:length(int)-1
val_ids=(gen>int(idx)&gen<=int(idx+1));
gen_log=[gen_log;mean(gen(val_ids))];
p_spk=[p_spk;mean(spks(val_ids))];
error_bar=[error_bar;sqrt(var(spks(val_ids)))/sqrt(sum(val_ids))];

end
idx=idx+1;
val_ids=(gen>int(idx));
gen_log=[gen_log;mean(gen(val_ids))];
p_spk=[p_spk;mean(spks(val_ids))];
error_bar=[error_bar;sqrt(var(spks(val_ids)))/sqrt(sum(val_ids))];

subplot(2,2,3);
errorbar(gen_log,p_spk,error_bar,'*');
hold on
x=[-1:0.01:1];
plot(x,N(x),'r');

%%

%% Generate response to stimulus



spkCondCollModel=struct('spksColl',[]);
spkCondCollModel(nConditions).spksColl=zeros(nTrials,size(condMovies{1},1));

for icond=1:nConditions
% Format movie into right format
mov_toformat=condMovies{icond};
mov=zeros(size(mov_toformat,2),size(mov_toformat,3),3,size(mov_toformat,1));

for itime=1:size(mov,4)
mov(:,:,1,itime)=squeeze(mov_toformat(itime,:,:));
mov(:,:,2,itime)=squeeze(mov_toformat(itime,:,:));
mov(:,:,3,itime)=squeeze(mov_toformat(itime,:,:));
end

% Generate response
ssta=ssta;
full_fit=full_fit;
s_use=full_fit(:,:,:,end:-1:1);
%s_use(:,:,:,15:end)=0; % DOUBT .. TODO ? Ask EJ
movie_time=size(mov,4);

gen_signals=zeros(movie_time,3);
for col=1:3;
    col
    st_temp=zeros(size(s_use,2),size(s_use,1),1,size(s_use,4)); % DOUBT .. Could be a reason for things to fail!!!!!
    for itime=1:30
        st_temp(:,:,:,itime)=s_use(:,:,col,itime)'; % DOUBT .. Could be a reason for things to fail!!!!!
    end
    s_use_new=st_temp;

Filtlen = size(s_use_new,4);
Filtdim1=size(s_use_new,1);
Filtdim2=size(s_use_new,2);

movie_new_len=movie_time;
mov2=zeros(Filtdim1 ,Filtdim2,movie_new_len+Filtlen-1);
mov2(:,:,Filtlen:movie_new_len+Filtlen-1)=squeeze(mov(:,:,col,:)); % Append zeros before the movie
sz=max(size(mov2,3)-size(s_use_new,4) + 1, 0);

gen_signals(:,col) = reshape(convn(mov2,squeeze(s_use_new(end:-1:1,end:-1:1,1,:)),'valid'),[sz,1]);
end
    %
gen=sum(gen_signals,2);
    
% gen is linear output
% N(gen) is after passing through non-linearity
pSpk = N(gen);



spksGen=cell(nTrials,1);
for iTrial=1:nTrials
spkCondCollModel(icond).spksColl(iTrial,:) = (poissrnd(pSpk)>0)';
end

spkCondCollModel(icond).spksColl=logical(spkCondCollModel(icond).spksColl);
end



figure;
for icond=1:nConditions
subplot(nConditions,1,icond);
plotSpikeRaster(spkCondCollModel(icond).spksColl,'PlotType','vertline');
title(sprintf('Data009 Vision ID: %d Avg Spk Rate: %f',InterestingCell_vis_id(ref_cell_number),0));
end

%%
nTrials1=nTrials;
figure;
for icond=1:nConditions
subplot(nConditions,1,icond);
[xPoints, yPoints]=plotSpikeRaster(spkCondColl(icond).spksColl,'PlotType','vertline');

plot(xPoints, yPoints+nTrials1, 'k');
spkCondColl(icond).xPoints=xPoints;
spkCondColl(icond).yPoints=yPoints;

hold on
[xPoints, yPoints]=plotSpikeRaster(spkCondCollModel(icond).spksColl,'PlotType','vertline');
plot(xPoints, yPoints, 'r');
spkCondCollModel(icond).xPoints=xPoints;
spkCondCollModel(icond).yPoints=yPoints;
title(sprintf('%s',cond_str{icond}));
end


figure;
for icond=1:nConditions
subplot(nConditions,1,icond);
xPoints = spkCondColl(icond).xPoints;
yPoints = spkCondColl(icond).yPoints;
nTrials1=max(yPoints(:));
plot(xPoints*120/20000, yPoints, 'k');

hold on
xPoints=spkCondCollModel(icond).xPoints;
yPoints=spkCondCollModel(icond).yPoints;
plot(xPoints, yPoints+max(yPoints(:)), 'r');
ylim([0,2*nTrials]);
title(sprintf('%s: data004 vis ID: %d ',cond_str{icond},InterestingCell_vis_id(ref_cell_number)));
end
