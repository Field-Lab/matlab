% Load data and STAs
startup_null_analyse_tenessee
%%


datafile = '2012-08-09-3/data002';
type_name= cell(1,1);
type_name{1}='On Parasol';

opt=struct('load_all',true);
datarun=load_data(datafile,opt)
datarun=load_sta(datarun)
datarun=load_params(datarun)
datarun=load_ei(datarun,'all','array_type',519);
%get_cell_ids(datarun,type_name) % Vision IDs - used for loading fitted
%STAs!

matlab_cell_ids=get_cell_indices(datarun,type_name);
stas=datarun.stas.stas(matlab_cell_ids);
n_cell=length(stas);

vision_id=1772; 
idx=[1:length(datarun.cell_ids)];
matlab_id=idx(datarun.cell_ids==vision_id);
cell_ei=datarun.ei.eis{matlab_id};
ssta=stas{matlab_cell_ids==matlab_id};


%%
addpath(genpath('~/Dropbox/Lab/Development/matlab-standard/code')); 
addpath(genpath('~/Dropbox/Lab/Development/matlab-standard/private/nishal/fwdfittingfunctions'));
fit_info = fit_sta(ssta)
%datarun = compute_sta_fits(datarun, 'all');

%%
%params=datarun.matlab.sta_fits{matlab_id}.initial_params;
params=fit_info.initial_params;
full_fit = sta_fit_function(params);

figure;
for itime=30:-1:1
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

    pause
    
end
%

%% Fit non Linearity
neuronPath = '/Volumes/Analysis/2012-08-09-3/data002/data002.neurons';
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuronPath);
spks=double(neuronFile.getSpikeTimes(vision_id));
TTLTimes= double(neuronFile.getTTLTimes());

%% Get Movie
mdf_file='/Volumes/Analysis/deprecated/movie-xml2/RGB-8-1-0.48-11111.xml'; % See in NOTEBOOKS
movie_time = 120*30*60;% TODO figure this out .. arbit no of frames written here! 
triggers=datarun.triggers;%[0:100/120:movie_time];
[mvi] = load_movie(mdf_file,triggers);
[~,height,width,duration,refresh] = get_movie_ath(mdf_file,...
    triggers, 1,2);
%%
var64=80;

mov=zeros(var64,40,3,movie_time);
delay=1;
    for itime=1:movie_time
        F=mvi.getFrame((itime-delay)/(120)).getBuffer;
        mov(:,:,1,itime)=double(reshape(F(1:3:end),width,height)); % Flip in each direction ? Doubt!
        mov(:,:,2,itime)=double(reshape(F(2:3:end),width,height)); % Flip in each direction ? Doubt!
        mov(:,:,3,itime)=double(reshape(F(3:3:end),width,height)); % Flip in each direction ? Doubt! % Check in stimulus rig and check 
        
    end
    
    mov=(mov-0.5)*(0.48/0.5);
%%

ssta=ssta;
full_fit=full_fit;
s_use=full_fit(:,:,:,end:-1:1);

gen_signals=zeros(movie_time,3);
for col=1:3;
    col
  st_temp=zeros(size(s_use,2),size(s_use,1),1,size(s_use,4)); % DOUBT .. Could be a reason for things to fail!!!!!
    for itime=1:30
        st_temp(:,:,:,itime)=s_use(:,:,col,end-itime+1)'; % DOUBT .. Could be a reason for things to fail!!!!!
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

    
    
% Fit non-linearity

%%
% Bin spikes in frame 
bin_spk=zeros(movie_time,1);
for ispk=1:length(spks)
numTriggers= sum(double(TTLTimes<=spks(ispk)));
EarlierTriggers = TTLTimes(TTLTimes<=spks(ispk));
lastTrig= EarlierTriggers(end);
frameNo=(numTriggers-1)*100 + (spks(ispk)-lastTrig)*120/20000;
bin_spk(floor(frameNo))=bin_spk(floor(frameNo))+1;
end

figure;
scatter(gen,bin_spk)
%%
    
addpath(genpath('../../gfield/'))
fit_params = fit_static_NL(bin_spk, gen)

N=@(x) exp(fit_params.b*x+fit_params.a);
inp=[-1:0.1:1];

hold on
plot(inp,N(inp),'r')
%% 
addpath(genpath('../../../code'));
datarun = get_snls(datarun, 'On Parasol','movie',mvi,'verbose',1,'marks','sigstix') % Try marks simple as well
%%
gen=datarun.stas.snls{matlab_id}.gen_signal;
spks=datarun.stas.snls{matlab_id}.spikes;

figure;
subplot(2,2,1);
hist(spks);
subplot(2,2,2);
scatter(gen,spks);
hold on;
x=[-1:0.01:1];
N=@(x) exp(datarun.stas.snls{matlab_id}.fit_params.a*x +datarun.stas.snls{matlab_id}.fit_params.b);
plot(x,N(x),'r');

subplot(2,2,4);
hist(gen,100);

int=[-1:0.1:1];
gen_log=[];
p_spk=[];
for idx=1:length(int)-1
val_ids=(gen>int(idx)&gen<=int(idx+1));
gen_log=[gen_log;mean(gen(val_ids))];
p_spk=[p_spk;mean(spks(val_ids))];
end
subplot(2,2,3);
plot(gen_log,p_spk);

%% 
fit_cell=cell(1,1);
fit_cell{1}=full_fit;%datarun.stas.stas{matlab_id};%full_fit;   

gen=compute_gen_signals(fit_cell,mvi);