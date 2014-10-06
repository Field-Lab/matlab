% Load data and STAs
startup_null_analyse_tenessee
%%


datafile = '2012-08-09-3/data000'%'2012-08-09-3/data002'; %'2012-08-21-1/data003';


opt=struct('load_all',true);
datarun=load_data(datafile,opt)
datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
datarun=load_params(datarun)
datarun = set_polarities(datarun);

%get_cell_ids(datarun,type_name) % Vision IDs - used for loading fitted
%STAs!

matlab_cell_ids=get_cell_indices(datarun,'all');
stas=datarun.stas.stas(matlab_cell_ids);
n_cell=length(stas);
% 
% vision_id=1772; 
% idx=[1:length(datarun.cell_ids)];
% matlab_id=idx(datarun.cell_ids==vision_id);
% cell_ei=datarun.ei.eis{matlab_id};
% ssta=stas{matlab_cell_ids==matlab_id};

%%
addpath(genpath('../../../code/'));
datarun = get_sta_summaries(datarun, datarun.cell_types{1}.name, ...
    'verbose',0,'keep_stas',0,'keep_rfs',1,'fig_or_axes',[],...
    'marks_params',struct( ...
        'strength','vector length', 'filter', fspecial('gauss',15,0.7), ...
        'thresh',5,'robust_std_method',1));
movie_spec='/Volumes/Analysis/movie-xml/BW-20-5-0.48-11111-30x30-60.35.xml'%' %RGB-8-1-0.48-11111.xml';%RGB-8-1-0.48-11111-80x40
datarun = load_java_movie(datarun, movie_spec);
datarun = get_snls(datarun, datarun.cell_ids(get_cell_indices(datarun, 'all')),'frames',-18:0,'stimuli',[],'new',true);
% TODO ask EJ .. 
%%
cellID=matlab_cell_ids(1);%matlab_id;
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
error_bar=[error_bar;sqrt(var(spks(val_ids)))];

for idx=1:length(int)-1
val_ids=(gen>int(idx)&gen<=int(idx+1));
gen_log=[gen_log;mean(gen(val_ids))];
p_spk=[p_spk;mean(spks(val_ids))];
error_bar=[error_bar;sqrt(var(spks(val_ids)))];

end
idx=idx+1;
val_ids=(gen>int(idx));
gen_log=[gen_log;mean(gen(val_ids))];
p_spk=[p_spk;mean(spks(val_ids))];
error_bar=[error_bar;sqrt(var(spks(val_ids)))];

subplot(2,2,3);
errorbar(gen_log,p_spk,error_bar,'*');
hold on
x=[-1:0.01:1];
plot(x,N(x),'r');

%% Better Load movie function?
mdf_file=movie_spec;%'/Volumes/Analysis/deprecated/movie-xml2/RGB-8-1-0.48-11111.xml';
triggers=datarun.triggers;
frames=20000%floor(max(triggers)*120) % 10 minutes
[mov,height,width,duration,refresh] = get_movie(mdf_file, triggers,frames);
mov=(mov-0.5);
%%
% Use Recursive Least Squares
% http://web.stanford.edu/class/ee263/lectures/06_ls-app.pdf - 6-21,
% Notation from there
delay=30;
mov_len=10000;
cellIds=matlab_cell_ids;
noCells=length(cellIds);

spk_coll = cell(noCells,1);

for icell=1:noCells
spk_coll{icell}=datarun.stas.snls{cellIds(icell)}.spikes;
end

filter_bank=cell(size(mov,1),size(mov,2),3);

for pixX=15%1:size(mov,1);
for pixY=12%1:size(mov,2);
for pixCol=1:1;
pixX
pixY
pixCol
mov_recons = squeeze(mov(pixX,pixY,pixCol,:));
q=zeros(delay*noCells+1,1);
P=zeros(delay*noCells+1,delay*noCells+1);
nSamples=mov_len-delay;
A=zeros(nSamples,delay*noCells+1);
b=zeros(nSamples,1);
for iFrame=1:nSamples%1800*120
    iFrame
    a=[1];% add a constant
    for icell=1:noCells
    a=[a;double(spk_coll{icell}(iFrame:1:iFrame+delay-1))];
    end
    a=full(a);
    y=mov_recons(iFrame);

    A(iFrame,:)=a';
    b(iFrame)=y;
    
    q=q+y*a;
    P=P+a*a';

end

 
filter=P\q;   
filter2=(A\(A'\(A'*b)));


% Make prediction
mov_pred=0*mov_recons;
recons_idx=mov_len:size(mov,4)-delay
for iFrame=recons_idx %1800*120
    iFrame
    a=[1];
    for icell=1:noCells
    a=[a;double(spk_coll{icell}(iFrame:1:iFrame+delay-1))];
    end
    a=full(a);
    
   mov_pred(iFrame)= filter'*a;
    
end

figure;
stairs(mov_recons(recons_idx(1)+100:recons_idx(1)+200),'b')
hold on
stairs(mov_pred(recons_idx(1)+100:recons_idx(1)+200),'r');
%hold on;
%stairs(double(mov_pred(recons_idx(1)+100:recons_idx(1)+200)>0)-0.5,'g');
xlim([1,100])
% correlation 

filter_bank{pixX,pixY,pixCol}.filter=filter;
filter_bank{pixX,pixY,pixCol}.correlation=corr(mov_pred(recons_idx),mov_recons(recons_idx));
corr(mov_pred(recons_idx),mov_recons(recons_idx))
corr(double(mov_pred(recons_idx)>0),mov_recons(recons_idx))
end
end
end



% save('/Volumes/Analysis/nishal/Reconstruct_2010-03-05-2_data018_1.mat','datarun');

% save('/Volumes/Analysis/nishal/Reconstruct_2010-03-05-2_data018_2.mat','filter_bank','mov','movie_spec','-v7.3');

%% 
corr_plot=zeros(26,size(mov,2));

iCol=1;
for xPix=1:size(mov,1)
    for yPix=1:size(mov,2)
    corr_plot(xPix,yPix)=filter_bank{xPix,yPix,iCol}.correlation;
    end
end
figure;
imagesc(corr_plot);
colormap gray

%% See Filter statistics

pixX=28;
pixY=12;
pixCol=1;

filters=zeros(noCells,delay);
filter_norm=zeros(noCells,1);
for icell=1:noCells
filters(icell,:)=filter_bank{pixX,pixY,pixCol}.filter(2+(icell-1)*delay:1+(icell)*delay);
filter_norm(icell)=norm(filters(icell,:));
end

figure;
subplot(2,1,1);
plot(filters');
subplot(2,1,2);
plot(filter_norm)

%% 
addpath(genpath('/Volumes/Analysis/nishal/cvx/'));
cvx_setup
%%
% Sparse filter
lambda = 10;
fLen=size(A,2);
cvx_begin
variable filter_sparse(fLen)
minimize ((0.5)*(A*filter_sparse-b)'*(A*filter_sparse-b) + lambda*norm(filter_sparse,1))
cvx_end

% Make prediction
mov_pred=0*mov_recons;
recons_idx=mov_len:108065-delay;
for iFrame=recons_idx %1800*120

    a=[1];
    for icell=1:noCells
    a=[a;double(spk_coll{icell}(iFrame:1:iFrame+delay-1))];
    end
    a=full(a);
    
   mov_pred(iFrame)= filter_sparse'*a;
    
end

figure;
stairs(mov_recons(recons_idx(1)+100:recons_idx(1)+200),'b')
hold on
stairs(mov_pred(recons_idx(1)+100:recons_idx(1)+200),'r');
%hold on;
%stairs(double(mov_pred(recons_idx(1)+100:recons_idx(1)+200)>0)-0.5,'g');
xlim([1,100])
% correlation 

filter_bank{pixX,pixY,pixCol}.filter=filter;
filter_bank{pixX,pixY,pixCol}.correlation=corr(mov_pred(recons_idx),mov_recons(recons_idx));
corr(mov_pred(recons_idx),mov_recons(recons_idx))
corr(double(mov_pred(recons_idx)>0),mov_recons(recons_idx))

%%
% Block Sparse filter
lambda = 10;
fLen=size(A,2);

cvx_begin
variable filter_sparse(fLen)

filters=cell(noCells,1);
for icell=1:noCells
filters{icell}=filter_sparse(2+(icell-1)*delay:1+(icell)*delay);
end

minimize ((0.5)*(A*filter_sparse-b)'*(A*filter_sparse-b) + lambda*norm(filter_sparse,1))
cvx_end

% Make prediction
mov_pred=0*mov_recons;
recons_idx=mov_len:108065-delay;
for iFrame=recons_idx %1800*120

    a=[1];
    for icell=1:noCells
    a=[a;double(spk_coll{icell}(iFrame:1:iFrame+delay-1))];
    end
    a=full(a);
    
   mov_pred(iFrame)= filter_sparse'*a;
    
end

figure;
stairs(mov_recons(recons_idx(1)+100:recons_idx(1)+200),'b')
hold on
stairs(mov_pred(recons_idx(1)+100:recons_idx(1)+200),'r');
%hold on;
%stairs(double(mov_pred(recons_idx(1)+100:recons_idx(1)+200)>0)-0.5,'g');
xlim([1,100])
% correlation 

filter_bank{pixX,pixY,pixCol}.filter=filter;
filter_bank{pixX,pixY,pixCol}.correlation=corr(mov_pred(recons_idx),mov_recons(recons_idx));
corr(mov_pred(recons_idx),mov_recons(recons_idx))
corr(double(mov_pred(recons_idx)>0),mov_recons(recons_idx))


