%% Train a bunch of neurons on real data
clear all;
close all;

nframes = 2^12;
filters_percenter = 3;
s = 13;
m = 40;

init_pars = struct('stim_type','lv',...
   'Mk',m,...
   'XsqK',0,...
   'wholesetsize', nframes,...
   'totalframes',72000,...
   'maxt',nframes,...
   'maxCGIter',80,...
   'maxNonsepIter',10,...
   'prelim_mode','nonsep',...
   'prelimIter',10,...
   'prelimFrames',nframes,...
   'datamode','real',...
   'preprocess_mode',@(x) x,...
   'fac',50,...
   'ps_fratio',20,...
   'cp_fratio',0,...
   'Nstep',@exp,'Nprime',@exp ,'Ndoubleprime',@exp,...
   'n',s^2,...
   'nofilters_ktime',m,...
   'ktime_basis',eye(m),...
   'k_spacing',pi/2,...
   'kspace_centers',[6 6; 6 7; 6 8; 7 6; 7 7 ; 7 8; 8 6; 8 7; 8 8],...
   'nofilters_kspace', 9*filters_percenter+1,...
   'min_std',1, 'max_std',4,....
   'nofilters_postspike',12,...
   'ps_spacing',pi,...
   'nofilters_coupling',0,...
   'bstretch',0.05, 'alpha_ps',-0.0045, 'alpha_cp',0,...
   'filtermode','sep_basis',...
   'hessmode','mult',...
   'padval',0.5,...
   'nocropflag',0,...
   'frame_offset',floor(200/(1/120)));

% monitor = [320 640]; <- unused
gopts = optimset('GradObj','on','Hessian','on','TolFun',10^(-8),'MaxIter',200,'display','iter');%,'DerivativeCheck','on','HessMult',@ll_hess_mult)

%% Dataset-specific stuff

% 2005-04-26 config
%load '~/glm_code/04-26-2005_cellids.mat';
%datafile = '~/data/shlens_newdata/2005-04-26-1-Spikes.mat';

% 2008-06-10 config
%load './06-10-2008_cellids.mat'
%datafile = '~/data/shlens_newdata/2008-06-10-0-Spikes.mat';

% 2008-08-27 config
%load '~/glm_code/08-27-2008_cellids.mat';
%datafile = '~/data/shlens_newdata/2008-08-27-6-Spikes.mat';

% Dec 2009 config
%datafile = '~/data/2009-12-03-2/2009-12-03-2-Spikes.mat';

% March 2010 config
datafile = '~/data/2010-03-05-1/spikes-2010-03-05-1.mat';
load(datafile)

switch(init_pars.stim_type)
   
   case 'wn'
      if (strfind(datafile,'2005'))
         stixel_val = 20;
      elseif (strfind(datafile,'2009-11-14'))
         stixel_val = 8;
      elseif(strfind(datafile,'2009-12-03'))
         stixel_val = 16;
      elseif(strfind(datafile,'2010-03-05'))
         stixel_val = 8;
      end
      init_pars.noise_type = 'bwn';
      X = get_wn_stimulus(datafile,init_pars.totalframes);
      
   case '1f'
      stixel_val = 10;
      init_pars.noise_type = 'g1fn';
      X = get_1f_stimulus(datafile);
      
   case 'mix'
      Xwn = get_wn_stimulus(datafile,floor(init_pars.totalframes/2));
      X1f = get_1f_stimulus(datafile);
      
      X = [Xwn X1f(:,1:(init_pars.totalframes - size(Xwn,2)))];
      
   case 'bl'
      stixel_val = 10;
      stim_params = struct('stixel', stixel_val, 'rf_size',8, 'rf_surround',4, ...
         'min_contrast', 0.01, 'lum_range', [0.5 0.5], ...
         'seed', 11111, 'integration_time', 0.25, ...
         'interval', 2);
   case 'lv'
      %load ~/lcvdata/shlens/data/repository/stim/Lava-2.0.mat;
      %load /users2/chaitu/lavastim_2009_12_3_72K.mat;
      %load /users2/chaitu/lavastim_2009_12_3_72K_resized.mat;
      %X = reshape(stim,160*80,size(stim,3));
      %stixel_val = 4;
      load /users2/chaitu/lava_v2_mu04_clip0025_binaryband_72000samples.mat;
      X = reshape(stim,80*160,72000); clear stim con band;
      stixel_val = 4;
end

load(datafile);
if 1
   switch(init_pars.stim_type)
      case '1f'
         spikes_lng = spikes_1f_lng;
         tstim = tstim_1f;
      case 'wn'
         %[blah blah lngidx] = intersect(cell_ids,cell_ids_wn_lng);
         %spikes_lng = subcellarray(spikes_wn_lng,lngidx);
         spikes_lng = spikes_wn_lng;
         tstim = tstim_wn_lng;
      case 'lv'
         %[blah blah lngidx] = intersect(cell_ids,cell_ids_lv_lng);
         %spikes_lng = subcellarray(spikes_lv_lng,lngidx);
         spikes_lng = spikes_lv_lng;
         tstim = tstim_lv_lng;
         
   end
end
init_pars.stim_height = 320/stixel_val;
init_pars.stim_width = 640/stixel_val;

% Resize the stimulus
%X = resize_stimulus(X,[80 160],[40 80],'bicubic');
%init_pars.n = 36; init_pars.stim_width = 80; init_pars.stim_height = 40;   init_pars.kspace_centers = [4 4; 3 3; 3 4; 4 3];


%% Single neurons
load(datafile)
%cell_idx = [265];
%cell_idx = 1:length(cell_ids);
%load ~/data/shlens_newdata/good_cell_idx.mat
%oad ~/data/2009-12-03-2/finalids.mat
%load ~/data/2009-12-03-2/adapt_ids.mat
%load ~/data/2009-12-03-2/best_ids.mat
%load ~/data/2009-12-03-2/ON_adapt_stable.mat
%load ~/data/2009-12-03-2/good_ids2.mat
%load ~/data/2009-12-03-2/ONOFF_goodids.mat
%cell_idx = find_idx(best_ids,cell_ids);
%cell_idx = find_idx(finalids,cell_ids);
%cell_idx = find_idx(ON_stable,c{ell_ids);
%cell_idx = find_idx(adapt_ids,cell_ids);
%cell_idx(1:2) = cell_idx([2 1]);
%load ~/data/2009-12-03-2/contrast_adaptation_stats.mat
%[sorted cell_idx] = sort(pinc_offset+pdec_onset,'descend');
%[blah blah cell_idx] = intersect(off_parasol,cell_ids);
%[blah blah cell_idx] = intersect(ON_goodids,cell_ids);
%offsets = floor([450 375 200 150 250]/tstim_lv_lng);
cell_idx = 1:length(cell_ids);

switch (init_pars.stim_type)
   case '1f'
      save_path = '/lcv/data/chaitu/glm/single_neuron/bestfits/1f_fits/';
      
      for j=1:length(cell_idx)
         [master_pstar basepars f g H] = multiple_neuron(cell_ids(cell_idx(j)),cell_ids(cell_idx(j)),init_pars,spikes_1f_lng,tstim_1f,cell_ids,init_pars.stim_type,[],save_path,gopts,X);
      end
      exit;
   case 'wn'
      save_path = '~/single_neuron/march10data/wn_fits/testing/';
      
      if (~isdir(save_path))
         mkdir(save_path); % Create the directory if it does not exist.
      end
      
      
      for j=148:length(cell_idx)
         init_pars.frame_offset = 400*120;
         init_pars = set_nframes(init_pars,2^15-40+1);
         leaveout_chunksize = 1;
         leaveout_nchunks = 30;
         leaveout_sepsize  = 1;
         init_pars.leaveout_idx = create_leaveout_idx(tstim,leaveout_chunksize,leaveout_sepsize,leaveout_nchunks,init_pars.maxt);
         %smooth_dt = 40;
         %init_pars.ext_timepts = 0:smooth_dt:(init_pars.maxt*tstim);
         [pstar basepars] = multiple_neuron(cell_ids(cell_idx(j)),cell_ids(cell_idx(j)),init_pars,spikes_lng,tstim,cell_ids,init_pars.stim_type,save_path,gopts,X);
      end
      exit;
   case 'lv'
      save_path = '~/single_neuron/march10data/lv_sb_fits/';
      if (~isdir(save_path))
         mkdir(save_path); % Create the directory if it does not exist.
      end
      
      
      for j=1:length(cell_idx)
         
         j
         
         init_pars.frame_offset = 40*120;
         init_pars = set_nframes(init_pars,1e4);
         leaveout_nchunks = 1;
         leaveout_chunksize = init_pars.maxt/2/120/leaveout_nchunks;
         leaveout_sepsize  = 2*leaveout_chunksize;
         init_pars.leaveout_idx = create_leaveout_idx(tstim,leaveout_chunksize,leaveout_sepsize,leaveout_nchunks,init_pars.maxt);
         init_pars.leaveout_idx = setdiff((1:init_pars.maxt)',init_pars.leaveout_idx);
         
         
         global CVLGP;
         global CVLGP_OLD;
         CVLGP = -1e8;
         CVLGP_OLD = -2e8;
         %init_pars.ext_timepts = 10;
         
         fprintf('The total range has %d frames. %d out of these are left out for cross-validation.\n',init_pars.maxt,length(init_pars.leaveout_idx));
         dx = zeros(init_pars.maxt,1);
         dx(init_pars.leaveout_idx) = true;
         plot(double(dx));
         % Make an external signal
         if 0
            %[x r] = plot_empirical_rate(spikes_lv_lng{find(cell_ids_lv_lng == cell_ids(cell_idx(j)))},tstim_lv_lng/50,tstim_lv_lng,600,3,'b');
            [x r] = plot_empirical_rate(spikes_lv_lng{find(cell_ids_lv_lng == cell_ids(cell_idx(j)))},tstim_lv_lng/50,1,600,10);
            init_pars.extrinsic = interp1(x,r,makeaxis(tstim_lv_lng,floor(x(end)/tstim_lv_lng)));
            plot(makeaxis(tstim_lv_lng,floor(x(end)/tstim_lv_lng)),init_pars.extrinsic,'r','LineWidth',2);
            init_pars.extrinsic = log(init_pars.extrinsic(:));
            init_pars.extrinsic = init_pars.extrinsic - mean(init_pars.extrinsic);
         end
         %hold on, plot_empirical_rate(spikes_lv_lng{find(cell_ids_lv_lng == cell_ids(cell_idx(j)))},tstim_lv_lng/50,0.1,600,1,'b');
         %pause;
         %close;
         %continue;
         
         if 0
            if (exist('S','var'))
               dataidx = find(S.cellids == cell_ids(cell_idx(j)));
               if (~isempty(dataidx))
                  load(S.filenames{dataidx});
                  init_pars.p0 = sep2nonsep(adjust_baserate_ext(master_pstar,basepars,stimpars,trainpars) + 0.01.*randn(size(master_pstar)),basepars);
                  
                  %init_pars.p0(get_pars_idx(basepars,1,1,'ps')) =  get_ortho_coeffs(init_pars.p0(get_pars_idx(basepars,1,1,'ps')),basepars.postspike_basis);
                  %kspaceidx = get_pars_idx(basepars,1,1,'kspace');
                  %ks1 = kspaceidx(1:basepars.nofilters_kspace); ks2 = kspaceidx(basepars.nofilters_kspace+1:end);
                  %init_pars.p0(ks1) = get_ortho_coeffs(init_pars.p0(ks1),basepars.kspace_basis);
                  %init_pars.p0(ks2) = get_ortho_coeffs(init_pars.p0(ks2),basepars.kspace_basis);
                  
                  clear basepars stimpars trainpars master_pstar;
                  init_pars.crop_idx = S.cpidces{dataidx};
               end
            end
         end
         
         % Set the nonlinearity
         %sp = spikes_lng{cell_idx};
         %[init_pars.Nstep init_pars.Nprime init_pars.Ndoubleprime NL{j}.coeffs NL{j}.range] = get_kx_NL(X,tstim,sp,nbins,2,init_pars.frame_offset+1:init_pars.frame_offset+init_pars.maxt);
         
         %init_pars = set_nframes(init_pars,8e3);
         tag = sprintf('ps_fratio%s_n%d_filtermode%s',num2str(init_pars.ps_fratio),sqrt(init_pars.n),init_pars.filtermode)
         [pstar basepars] = multiple_neuron(cell_ids(cell_idx(j)),cell_ids(cell_idx(j)),init_pars,spikes_lng,tstim,cell_ids,init_pars.stim_type,save_path,gopts,X,tag);
         
      end
      
      
      % Evaluate the fits
      evaluate_singlecells(save_path,[],[],[],[],X,[]);
      
      
      % Run varswitch
      varswitch_eval3(save_path,[],[],[],1000,[]);
      
      %exit;
      
      
      
end

%% Fit with fixed filter

S = recover_pars('~/single_neuron/dec03data/best_ids/lv_nsb_fits/');
%load ~/data/2009-12-03-2/adapt_ids.mat;
load ~/data/2009-12-03-2/best_ids.mat
cell_idx = find_idx(best_ids,cell_ids);

save_path = '~/single_neuron/dec03data/best_ids/lv_nsb_reopt_fits/';

gopts = optimset('GradObj','on','Hessian', 'on', 'TolFun',10^(-6),'MaxIter',25,'display','iter');%,'DerivativeCheck','on','HessMult',@ll_hess_mult)

for j=1:length(cell_idx)
   
   % Just reopt
   load(S.filenames{find(S.cellids == cell_ids(cell_idx(j)))});
   basepars.p0 = master_pstar;
   basepars.crop_idx = basepars.crop_idx{1};
   [pstar basepars] = multiple_neuron(cell_ids(cell_idx(j)),cell_ids(cell_idx(j)),basepars,spikes_lng,tstim,cell_ids,init_pars.stim_type,save_path,gopts,X);
   
   
   % Reopt with fixed filt
   %kidx = find(S.cellids == cell_ids(cell_idx(j)));
   % Obtain the filter
   %init_pars.K{1} = S.K{kidx};
   % Obtain crop indices
   %init_pars.crop_idx = S.cpidces{kidx};
   %init_pars.frame_offset = ON_stable_offsets(j);
   %[pstar basepars] = multiple_neuron(cell_ids(cell_idx(j)),cell_ids(cell_idx(j)),init_pars,spikes_lng,tstim,cell_ids,init_pars.stim_type,save_path,gopts,X);
   
   
end

evaluate_singlecells(save_path,[],[],[],X);
clear X;
varswitch_eval3(save_path,[],[],[]);
exit;





%% Further optimization of nonlinearity after fitting parameters

oldparspath = '~/single_neuron/dec03data/saids/lv_nsb_ext_fits/';
S = recover_pars(oldparspath);
save_path = '~/single_neuron/dec03data/saids/lv_nsb_ext_nl_fits/';

load ~/data/2009-12-03-2/strongadapt_ids.mat;
cell_idx = find_idx(strongadapt_ids,cell_ids);

NITER = 1;% number of iterations
nl_nbins = 400; % number of bins used to fit nonlinearity
nl_deg = 2;
master_NL = cell(length(S.cellids),1);
%master_P = cell(length(S.cellids),1);

% Setup options
gopts = optimset('GradObj','on','Hessian', 'on', 'TolFun',10^(-4),'MaxIter',20,'display','iter');%,'DerivativeCheck','on','HessMult',@ll_hess_mult)

for j=1:length(S.cellids)
   
   % Reset anything that might have changed
   init_pars.Nstep = @exp; init_pars.Nprime = @exp; init_pars.Ndoubleprime = @exp;
   if (isfield(init_pars,'p0'))
      init_pars = rmfield(init_pars,'p0');
   end
   if (isfield(init_pars,'crop_idx'))
      init_pars = rmfield(init_pars,'crop_idx');
   end
   
   
   % Set the frame offset and range
   %best_idx = find(best_ids == S.cellids(j));
   %init_pars.frame_offset = best_offsets(best_idx);
   %init_pars = set_nframes(init_pars,best_range(best_idx));
   
   % Load the initial fit
   load(S.filenames{j});
   pstar = master_pstar;
   init_pars = basepars;
   
   % Run the initial fit with exponential nonlinearity
   %[pstar basepars] = multiple_neuron(cell_ids(cell_idx(j)),cell_ids(cell_idx(j)),init_pars,spikes_lng,tstim,cell_ids,init_pars.stim_type,'./',gopts,X,'tmp.mat');
   
   if 0
      % Convert to fix filt mode
      init_pars.filtermode = 'fixfilt';
      init_pars.K{1} = S.K{j};
      pstar = [pstar(get_pars_idx(basepars,1,1,'ext'));...
         pstar(get_pars_idx(basepars,1,1,'b'));...
         1;... % Knorm
         pstar(get_pars_idx(basepars,1,1,'ps'));...
         pstar(get_pars_idx(basepars,1,1,'cp'));];
   else
      % Stay in sep_basis mode, no preprocessing needed
   end
   
   
   
   % Iterate between nonlinearity and parameter estimation
   NL = zeros(nl_deg+1,NITER);
   sfile = S.filenames{j};
   
   sp{1} = spikes_lng{find(cell_ids == S.cellids(j))};
   
   for k=1:NITER
      
      % Re-estimate nonlinearity
      [fnew f master_NL{j}.coeffs master_NL{j}.range] = sub_NL_ll_eval_wrapper(sfile,X,sp{1},nl_nbins,nl_deg);
      close all;
      
      if 0%(((fnew-f)/f*100 < 0.1) && k > 1) % If the nonlinearity adjustment does not improve LL, dont refit. quit.
         fprintf('Nonlinearity optimization has not increased log-likelihood of training data by more than 1 percent. Quitting.\n');
         break;
      end
      
      % Reset the nonlinearity - make it piecewise on the range, constant elsewhere
      init_pars = set_nl(init_pars,master_NL{j}.coeffs,master_NL{j}.range);
      
      if (k== NITER)
         spath = save_path; % Save only on last iteration
         sfile = [];
      else
         spath = './';
         sfile = 'tmp.mat';
      end
      
      % Set the initialization to the last fit
      init_pars.p0 = pstar;
      init_pars.crop_idx = basepars.crop_idx{1};
      
      % Run the fitting
      [pstar basepars] = multiple_neuron(S.cellids(j),S.cellids(j),init_pars,spikes_lng,tstim,cell_ids,init_pars.stim_type,spath,gopts,X,sfile);
      %P(:,k) = pstar;
      
      close;
   end
   
   
end

save(strcat(make_outpath(save_path),'nlparsinfo.mat'),'master_NL','NITER','nl_nbins','nl_deg','-v7.3');

% Run evaluations
evaluate_singlecells(save_path,[],struct('nonlinflag',1),datafile,X);
clear X;
varswitch_eval3(save_path,[],[],[],500,1);
exit;

%% Single neurons - bootstrap


load(datafile);
switch(init_pars.stim_type)
   case '1f'
      save_path = '/lcv/data/chaitu/glm/single_neuron/newfits/1f_nonsep_bootstrap/ON/';
      spikes_lng = spikes_1f_lng;
      tstim = tstim_1f;
   case 'wn'
      save_path = '/lcv/data/chaitu/glm/single_neuron/newfits/wn_nonsep_bootstrap/ON/';
      spikes_lng = spikes_wn_lng;
      tstim = tstim_wn;
end

nbootstraps = 25;


[blah blah onidx] = intersect(on_parasol,cell_ids(cell_idx));
[blah blah offidx] = intersect(off_parasol,cell_ids(cell_idx));

cell_idx = cell_idx(onidx);


cells_totrain = cell_ids(cell_idx);

results = cell(length(cells_totrain),1);
base_blocksize = ceil(((2^13)-init_pars.Mk+1));
for j=1:length(cells_totrain)
   
   switch(init_pars.stim_type)
      case 'wn'
         blocksize = base_blocksize;
      case '1f'
         blocksize = floor(base_blocksize*length(spikes_wn_lng{cell_idx(j)})/length(spikes_1f_lng{cell_idx(j)})) % adjust to get approx. same number of spikes
   end
   
   nblocks = floor((init_pars.totalframes - init_pars.frame_offset)/blocksize)
   results_j = bootstrap_neuron(init_pars,nbootstraps,cells_totrain(j),nblocks,blocksize,X,'',init_pars.stim_type,cell_ids,spikes_lng,tstim,gopts);
   results{j} = results_j;
   id = cells_totrain(j);
   if (~isempty(results_j))
      save(strcat(save_path,sprintf('bootstrap_summary_cell_%d.mat',cells_totrain(j))),'results_j','init_pars','nblocks','blocksize','tstim','nbootstraps','save_path','datafile','id');
   end
end

save(strcat(save_path,'bootstrap_summary_allcells.mat'),'results','init_pars','nblocks','blocksize','tstim','cells_totrain','save_path','datafile','nbootstraps');


%% 2009 Data set - single neurons

save_path = '/lcv/data/chaitu/glm/single_neuron/2009_dataset/blrep_fits/';

for j=1:length(cell_idx)
   multiple_neuron(cell_ids(cell_idx(j)),cell_ids(cell_idx(j)),init_pars,spikes_lng,tstim_bl,cell_ids,init_pars.stim_type,stim_params,save_path,gopts,Xlong);
end


%% 2001 Data set - single neurons - separate stimulus

% Set up the stimulus/spikes
load ~/2001data/2001-03-16-1-spikes.mat;

ncells = length(spikes04);
gopts = optimset(gopts,'MaxIter',50);
cell_ids = 1:ncells;
init_pars.stim_width =1 ;
init_pars.stim_height = 1;
init_pars.n = 1;

for j=1:4
   
   X = evalin('base',strcat('X0',num2str(j+3)));
   spikes_j = evalin('base',strcat('spikes0',num2str(j+3)));
   save_path = strcat('/lcv/data/chaitu/glm/single_neuron/2001_dataset/wn0',num2str(j+3),'_fits/')
   
   for k=1:ncells
      
      nframes = ceil(spikes_j{k}(end)/tstim_wn)
      
      init_pars.wholesetsize = nframes;
      init_pars.totalframes = nframes;
      init_pars.maxt = nframes;
      init_pars.prelimFrames = nframes
      init_pars.prelimIter = 0;
      
      multiple_neuron(cell_ids(j),cell_ids(j),init_pars,spikes_j,tstim_wn,cell_ids,init_pars.stim_type,wn_stim_params,save_path,gopts,X);
      
   end
   
   
end


%% 2001 Data set - single neurons - mixed stimulus

% Set up the stimulus/spikes
load ~/2001data/2001-03-16-1-spikes.mat;

ncells = length(spikes04);


frames_per_condition = 50000;
nconditions = 4;
X = zeros(size(X04,1),nconditions*frames_per_condition);
time_per_condition = frames_per_condition * tstim_wn;

init_pars.wholesetsize = frames_per_condition*nconditions;
init_pars.totalframes = frames_per_condition*nconditions;
init_pars.maxt = frames_per_condition*nconditions;
init_pars.prelimFrames = 1;%frames_per_condition*nconditions;
init_pars.prelimIter = 0;


%init_pars.rms_mlength = 0;
%init_pars.rms_slength = 100;

save_path = '/lcv/data/chaitu/glm/single_neuron/2001_dataset/wnmix_fits/50K_frames_softrect_fsf01/';
spikes = cell(ncells,1);
for k=1:ncells
   spikes{k} = [];
end

for j=1:nconditions
   Xtmp = evalin('base',strcat('X0',num2str(j+3)));
   X(:,(j-1)*frames_per_condition+1:j*frames_per_condition) = Xtmp(:,1:frames_per_condition);
   
   spikes_j = evalin('base',strcat('spikes0',num2str(j+3)));
   
   for k=1:ncells
      spikes{k} = [spikes{k}; ((j-1)*time_per_condition+spikes_j{k}(spikes_j{k} < time_per_condition))];
   end
   
end

cell_ids = 1:6;
init_pars.stim_width =1 ;
init_pars.stim_height = 1;
init_pars.n = 1;
init_pars.frame_offset = 0;
gopts = optimset(gopts,'MaxIter',50);


for j=1:length(cell_ids)
   multiple_neuron(cell_ids(j),cell_ids(j),init_pars,spikes,tstim_wn,cell_ids,init_pars.stim_type,wn_stim_params,save_path,gopts,X);
end



%% Multiple neurons  - Create multiple cell index matrix

load(datafile)
%load /lcv/data/chaitu/1f_72Kframes_08_2008.mat;
%cell_list = off_parasol;
cell_list = cell_ids;
load ~/data/2009-12-03-2/adapt_ids.mat
load ~/data/2009-12-03-2/dist_table.mat
dist_table(abs(dist_table) < eps) = inf;
cell_table = dist_table < 6;
%cell_table = table_off_parasol;

ncell = 8;

multicell_ids = zeros(size(cell_list,1),ncell);

%good_idx = logical(true(size(cell_list,1),1));

for j=1:size(multicell_ids,1)
   
   %if (find(multicell_ids(:) == cell_list(j)))
   %    continue;
   %end
   %good_idx(j) = 1;
   
   % Retrieve j'th cell
   adj = find(cell_table(j,:));
   num = min(ncell-1,length(adj));
   multicell_ids(j,1) = cell_list(j);
   multicell_ids(j,2:num+1) = cell_list(adj(1:num));
   
   cellrow = multicell_ids(j,:)';
   
   if (find(multicell_ids(1:j-1,:) == cellrow(1)))
      good_idx(j) = false;
   end
   
   
end

%multicell_ids = multicell_ids(good_idx,:)

%multicell_ids = multicell_ids(floor(size(multicell_ids,1)/2):end,1)

%multiple_neuron(cell_ids(cell_idx),cell_ids(cell_idx),init_pars,datafile,init_pars.stim_type,stim_params,1,gopts);
[adapt_on_ids blah blah] = intersect(on_parasol,adapt_ids);
[blah blah cell_idx] = intersect(adapt_on_ids, cell_ids);
multicell_ids = multicell_ids(cell_idx,:);


save_path = '~/multiple_neuron/dec03data/adapt_ids/lv_newsepbasis_fits/';

for j=1:size(multicell_ids,1)
   
   cellrow = multicell_ids(j,:);
   lastidx = find(cellrow == 0,1);
   if (isempty(lastidx))
      lastidx = size(multicell_ids,2)+1;
   end
   
   cellrow = cellrow(1:lastidx-1)
   
   multiple_neuron(cellrow(1),cellrow(:),init_pars,spikes_lng,tstim,cell_ids,init_pars.stim_type,[],save_path,gopts,X);
end

exit;


%% Multiple neurons - create simulated data

pars_type = 'wn';
init_pars.stim_params = stim_params;
load(datafile);

% Load the parameters
switch (pars_type)
   case 'wn'
      load '/lcv/data/chaitu/glm/single_neuron/wn_04_29_2009_12by12_good/cell181_alldata_12by12_wn_nonsep_27-Apr-2009.mat';
   case '1/f'
      load '/lcv/data/chaitu/glm/single_neuron/1f_04_29_2009_12by12_good/cell181_alldata_12by12_1f_nonsep_27-Apr-2009.mat';
end
multicell_ids = 181;

% Load the stimulus
switch (init_pars.stim_type)
   case 'wn'
      X = export_stimulus(init_pars.stim_type, init_pars.totalframes, [], init_pars.stim_params)';
      trainpars.dt = tstim_wn/basepars.fac;
   case '1/f'
      load /lcv/data/chaitu/1f.mat;
      trainpars.dt = tstim_1f/basepars.fac;
end

switch (init_pars.preprocess_mode)
   case 'log'
      X = log(X+eps);
   case 'raw'
      % do nothing
end


% Adjust parameters if the stimulus does not match conditions under which
% parameters were trained

if (strcmp(pars_type,'wn') && strcmp(init_pars.stim_type,'1/f'))
   % Make the resolution of the parameters finer
   ptheta = reshape(repcols(reshape(pstar(2:1+basepars.n*basepars.Mk,1),basepars.n,basepars.Mk),init_pars.Mk/basepars.Mk),init_pars.Mk*basepars.n,1);
   pstar = [pstar(1);ptheta;pstar(end-basepars.nofilters_postspike+1:end,1)];
end

if (strcmp(pars_type,'1/f') && strcmp(init_pars.stim_type,'wn'))
   % Make the resolution of the parameters coarser
   ptheta = reshape(pstar(2:1+basepars.n*basepars.Mk,1),basepars.n,basepars.Mk);
   ptheta = ptheta(:,1:basepars.Mk/init_pars.Mk:end); %(take every other frame)
   ptheta = ptheta(:);
   pstar = [pstar(1);ptheta;pstar(end-basepars.nofilters_postspike+1:end,1)];
end

% Copy relevant parameters
basepars.Mk = init_pars.Mk;
basepars.maxt = 2^14-basepars.Mk+1;
basepars.datamode = 'simulate';
basepars.preprocess_mode = init_pars.preprocess_mode;
basepars.prelim_mode = 'sep';
basepars.noise_type = init_pars.noise_type;
basepars.wholesetsize = basepars.maxt;
basepars.totalframes = init_pars.totalframes;
basepars.prelimIter = 5;
basepars.prelimFrames = basepars.maxt;
basepars.stim_type = init_pars.stim_type;

[spikes lcifs] = simulate_model(basepars,pstar(:,1),X(:,1:basepars.maxt),[],trainpars.dt,1);
sim_data{1} = find(spikes)*trainpars.dt;



multiple_neuron(multicell_ids(:,1),multicell_ids,basepars,datafile,init_pars.stim_type,init_pars.stim_params,1,gopts,sim_data,pstar(:,1));


