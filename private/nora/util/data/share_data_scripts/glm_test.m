piece = '2013-08-19-3';
stim_type = 'NSEM';
cell = 5;
tstim = .00832750;

%% loading and organizing stimulus
load(['/Volumes/Lab/Users/Nora/data_files/Data/CarlosData/' stim_type '-' piece '-CellData.mat'])
load(['/Volumes/Lab/Users/Nora/data_files/Data/CarlosData/' stim_type '-' piece '-StimData.mat'])
eval(['blockedmoviecell = ' stim_type 'StimData.FitMovie;'])
eval(['testmovie = ' stim_type 'StimData.TestMovie;'])
eval(['cells = fieldnames(' stim_type 'CellData);'])
clear WNStimData
% concatenate the fit movie (different blocks
n_blocks = length(blockedmoviecell);
height       = size(blockedmoviecell{1},2);
width        = size(blockedmoviecell{1},1);
fitframes    = size(blockedmoviecell{1},3);
totalframes       = n_blocks * ( fitframes) ;
concat_fullfitMovie = uint8(zeros(width, height, totalframes)) ;
for i_blk = 1:n_blocks
        framenums = ( (i_blk -1)*fitframes + 1 ) :  (i_blk *fitframes);  
        concat_fullfitMovie(:,:,framenums) = blockedmoviecell{i_blk};    
end
fitmovie = concat_fullfitMovie;
clear concat_fullfitMovie blockedmoviecell framenums height i_blk totalframes width

%% cell specific
eval(['names = fieldnames(' stim_type 'CellData);'])
eval(['spikes = ' stim_type 'CellData.' names{cell} '.Spikes;'])
eval(['STA = ' stim_type 'CellData.' names{cell} '.STA;'])

%% concatenate spikes
testspikes = spikes(1:2:end);
blockedspikes = spikes(2:2:end);
t_start   = 0;
T_SP = []; blk_count = 0;
dur = tstim * fitframes;
for k = 1:n_blocks
	blk_count = blk_count + 1;
	t_sp_full = blockedspikes{k} ; % unit of time: sec, 0 for the onset of the block
	t_sp      = t_sp_full(find(t_sp_full >  t_start));
	t_sp = t_sp - t_start;
	t_spcontext = t_sp + ( blk_count -1 )*dur;
	T_SP = [T_SP ; t_spcontext];
end
fitspikes = T_SP;
clear T_SP blk_count blockedspikes dur fitframes k n_blocks spikes t_sp t_sp_full t_spcontext t_start 

%%
[~,center] = STA_Test(fitspikes, fitmovie, 1, tstim);
fittedGLM = glm_fit(fitspikes, fitmovie, center, 'WN_STA', STA, 'monitor_refresh', 1/tstim);
fittedGLM.xval = glm_predict(fittedGLM,testmovie, 'testspikes', testspikes);
figure(1); plotfilters(fittedGLM)
figure(2); plotrasters(fittedGLM.xval, fittedGLM)
