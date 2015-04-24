% Function for fitting multiple neurons

% Inputs: 

% cells_totrain - list of cell ids to fit
% cells_touse - list of cell ids whose spikes to use
% basepars - struct with basic paramters (see fitting_script.m)
% datafile - file where spike data is located
% stim_type - string specifying stimulus type "wn" or "1/f" or "bl"
% stim_params - stimulus parameters for stimulus generation (white nosie)
% save_dir - path of directory in which to save output
% gopts - any optimization options you want to pass on to fminunc

function [master_pstar basepars master_f master_g master_H]  = multiple_neuron_c(cells_totrain,cells_touse,basepars,spikes_lng,tstim,cell_ids,stim_type,save_dir,gopts,X,C,tag)

% Set the filter mode ('sep' is rank-2, 'nonsep' is full rank)
if (~isfield(basepars,'filtermode'))
    filtermode = 'sep_raw'
else
    filtermode = basepars.filtermode
end

stimpars.dt = tstim; % this is the time (in sec) each stimulus frame is presented for

if (~isfield(basepars,'preprocess_mode'))
    basepars.preprocess_mode = @(x) x;
end

if (~exist('C','var'))
    C = [];
end

X = basepars.preprocess_mode(X); % in case there are 0s add an epsilon to X


% Form the file suffix for saving output
file_suffix = strcat(stim_type,'_',filtermode,'_','_',date);

% Base parameters (used by training procedure)
[sdim tdim ddim] = size(X); % space/time dimension of stimulus
basepars.stim_n = sdim;

multiflag = (ddim>1);

[basepars trainpars]  = create_bases_fns(basepars,stimpars,trainpars);
    
% Automatically figure out sampling times if needed
if (isfield(basepars,'ext_timepts') && ~isempty(basepars.ext_timepts))

    if (strcmp(basepars.ext_timepts,'auto'))
        % Find the cutoff frequency of the postspike basis in Hertz
        fcutoff = 1/(basepars.ps_fratio*basepars.Mk*stimpars.dt);
        extsample_rate = (1/fcutoff)/stimpars.dt; % # of stimulus bins per sample    
    else
        extsample_rate  = floor(basepars.ext_timepts(1)/stimpars.dt);
    end
    
    fprintf('Modeling extrinsic signal at a sampling rate of %0.3f sec per sample.\n',extsample_rate*stimpars.dt);
    
    basepars.ext_timepts = (1:extsample_rate:basepars.maxt)';
    nsamples = length(basepars.ext_timepts);

    % Take samples
    samples = randn(nsamples,1);% x(sample_idx);

    % Interpolate
    [y b] = interp(samples,extsample_rate,4,0.15);

    % Save the interpolating filter
    basepars.interpfilt = b;
    
end

% Set the number of neurons TO TRAIN and those TO USE (Neff >= N always)
N = basepars.Nneurons;
Neff = length(cells_touse); % 'effective number of neurons'

fprintf('N=%d\tNeff=%d\n',N,Neff);
sp_times = cell(Neff,1);

% Obtain the spike times for each cell that is being used in training and
% store in a cell array as well as in a binary matrix
spikeT = basepars.fac*basepars.maxt;

if (~isfield(basepars,'frame_offset'))
    stim_offset = 0;
else
    stim_offset = stimpars.dt*basepars.frame_offset;
end

if (~multiflag)
    
    trainpars.D = sparse(logical(false(spikeT,Neff)));
    duration = stimpars.dt * basepars.maxt; % amount of time in seconds of stimulus segment for training
    
    [sp_times trainpars.D trainpars.negSpikes] = pull_spikes(spikes_lng,cell_ids,cells_touse,spikeT,duration,trainpars.dt,stim_offset);
    
    % Plot the training segment
    %figure(4), plot_empirical_rate(sp_times,trainpars.dt,0.1,duration,0.1,'b');
    
    1;
    
    % Compute the postspike/coupling basis gradients (convolve each spike train with the basis fns)
    fprintf('Computing the spike train convolutions...');
    spikeT_offset = basepars.frame_offset*basepars.fac;
    trainpars.psbasisGrad = grad_basis([trainpars.negSpikes; trainpars.D],basepars.postspike_basis,spikeT_offset);
    if (Neff > 1)
        trainpars.cpbasisGrad = grad_basis([trainpars.negSpikes; trainpars.D],basepars.coupling_basis,spikeT_offset);
    else
        trainpars.cpbasisGrad = [];
    end
    fprintf('done.\n');
else
    sp_times = pull_spikes(spikes_lng,cell_ids,cells_touse,0,600,trainpars.dt,stim_offset);
    trainpars.D = ones(1,Neff);
    trainpars.analog_spikes = sp_times;
end

% Initialize master variables

% Set the number of parameters per neuron based on the filter mode
npars_perneuron = get_npars(basepars,Neff);

master_pstar = [];
master_p0 = [];
master_f = zeros(N,1);
master_g = [];


% Set the stimulus
frame_idx = basepars.frame_offset+1:basepars.frame_offset+basepars.maxt;
if(multiflag)
    for j=1:size(X,3)
        stimpars.x(:,:,j) = X(:,frame_idx,j);
    end
else
    stimpars.x = X(:,frame_idx);
end
%clear X;

1;

saveFlag = ~isempty(save_dir) && ~strcmp(save_dir,'');

% Make the output directory
if (saveFlag)
    cell_prefix = strcat('cells_',form_cell_prefix(cells_totrain(:)));
    if (N > 1)
        fprintf('Creating a folder %s to save output\n',cell_prefix);
        mkdir(save_dir,cell_prefix);
        save_dir = strcat(save_dir,cell_prefix,'/');
    end
end

1;

% Iterate through each cell and train

master_crop_idces = cell(N,1);%zeros(basepars.n,N);

for j=1:N

    % Reorder the spike matrix so this neuron's col comes first (convention);
    %neworder = [j 1:j-1 j+1:Neff];
    %trainpars.D = D(:,neworder);
    %trainpars.psbasisGrad = permute_cellarray(PSBCONV,neworder);
    %trainpars.cpbasisGrad = permute_cellarray(CPBCONV,neworder);
    
    % Set the number of neurons to be trained at this time to 1 (fit one at a time insteady of simultaneously)
    basepars.Nneurons = 1;
    trainpars.baseneuron_idx = [j]; % index of the neuron(s) being trained simultaneously
    
    1;
    % Get the STA and initialization point
    
    if (~strcmp(basepars.filtermode,'fixfilt') && ~isfield(basepars,'p0'))

        [stas rstas cpidces p0 box_m box_n] =  find_init_pt(sp_times{trainpars.baseneuron_idx},basepars,stimpars,trainpars,basepars.prelimIter,basepars.prelimFrames,0);
        
        1;
        
        if (isempty(cpidces))
            fprintf('Cell id %d has an invalid cropping window!\n',cells_totrain(j));
            return; % Skip this cell!
        end
    elseif (strcmp(basepars.filtermode,'fixfilt') && ~isfield(basepars,'p0'))
        p0 = 0.01.*randn(npars_perneuron*basepars.Nneurons,1);
    else
        fprintf('Using supplied initialization\n');
        p0 = basepars.p0;
    end

    
    if (~isfield(basepars,'crop_idx'))
        if (size(cpidces,1) ~= basepars.n)
            fprintf('multiple_neuron: basepars.n has been changed to %d\n',size(cpidces,1));
            if (strcmp(basepars.filtermode,'sep_basis'))
                fprintf('ERROR: kspace-basis is not ready to handle truncated RFs!Choose another filtermode.\n');
                return;
            end
            basepars.n = size(cpidces,1);
        end
        basepars.box_m = box_m;
        basepars.box_n = box_n;
        
        %basepars.ACmat = ACmat;
        basepars.crop_idx = cpidces;
    else
        fprintf('Using supplied cropping indices\n');
    end
        
    
    master_crop_idces{j} = basepars.crop_idx;

    timeoffset = basepars.frame_offset*stimpars.dt; % offset time in seconds 
    fprintf('Using frames %d-%d spanning the time range %f-%f sec.\n',basepars.frame_offset+1,basepars.frame_offset+basepars.maxt,timeoffset,timeoffset+duration);
    
    fprintf('Training cell %d/%d with id %d. The neuron is coupled to cells :',j,length(cells_totrain),cells_totrain(j));
    for r=[1:j-1 j+1:length(cells_touse)]
        fprintf('%d, ',cells_touse(r));
    end
    fprintf(' Npar: %d\n',npars_perneuron);
    
    1;
    
    switch filtermode % initialize the rest of the params based on the filter mode
        case 'sep_basis'
        case 'sep_raw'
        case 'nonsep'
            trainpars.lgrad = nonsep_lgrad(basepars,stimpars,trainpars);
        case 'fixfilt'
    end
    
    % Reset any frozen parameters
    if (isfield(basepars,'frozen_idx') && isfield(basepars,'frozen_vals') && ~isempty(basepars.frozen_idx))
        fprintf('Freezing parameters at %d indices...\n',length(basepars.frozen_idx));
        p0(basepars.frozen_idx) = basepars.frozen_vals;
    end
    
    % Start the training process
    1;
    
    if (isfield(basepars,'stimreg') && strcmp(basepars.stimreg,'L2') && isfield(basepars,'lambda_reg') && basepars.lambda_reg > 0)
        fprintf('L2-Regularizing stimulus filter with lambda=%0.5f\n',basepars.lambda_reg);
    end
    
    
1;
    
    switch filtermode
        case 'sep_basis'
            [pstar fstar gstar] = train_glm2(p0,basepars,stimpars,trainpars,gopts,1);
        case 'sep_raw'
            [pstar fstar gstar] = train_glm2(p0,basepars,stimpars,trainpars,gopts,1);
        case 'nonsep'
            [pstar fstar gstar] = train_glm_nonsep(p0,basepars,stimpars,trainpars,gopts);
        case 'fixfilt'
            [pstar fstar gstar] = train_glm_fixfilt(p0,basepars,stimpars,trainpars,gopts,1);
        case 'nonsep_cg'
            [pstar fstar gstar] = train_glm_nonsep_cg(p0,basepars,stimpars,trainpars,gopts);
    end
   
    
    % Reset any frozen parameters
    if (isfield(basepars,'frozen_idx') && isfield(basepars,'frozen_vals') && ~isempty(basepars.frozen_idx))
        pstar(basepars.frozen_idx) = basepars.frozen_vals;
    end
    
    % Concatenate the results
    1;
    %param_idx = (j-1)*npars_perneuron+1:j*npars_perneuron;
    master_p0 = [master_p0; p0];%(param_idx) = p0;
    master_pstar = [master_pstar; pstar];%(param_idx) = pstar;
    master_f(j) = fstar;
    master_g = [master_g; gstar];%(param_idx) = gstar;
    if (nargout >= 5)
        [blah1 blah2 master_H] = ll_func2(pstar,basepars,stimpars,trainpars);
    end
    
    
    1;
    
    if 0;%(saveFlag)
        close all;
        
        figure(2);
        h = figure(2);
        set(h,'position',[0 0 1920 1072]);
        
        % Set the name of the file to save output
        figname = strcat(save_dir,'cell_',num2str(cells_totrain(j)),'_cpl_',form_cell_prefix(cells_touse([1:j-1 j+1:Neff])),'_',file_suffix);        
        1;
        plot_opt_results(pstar,p0,basepars,stimpars,trainpars,filtermode,fstar,gstar);
        
        savefig(figname,h,'eps');
        %print(h,'-dpdf',strcat(fileprefix,'cell',num2str(cells_totrain(j)),'_',filesuffix));
        
               
        close all;
        
    end
end



if (saveFlag)

        % Remove any space-taking variables before saving
    
        if (isfield(stimpars,'x'))
            stimpars = rmfield(stimpars,'x');
        end
        
        if (isfield(trainpars,'psbasisGrad'))
            trainpars = rmfield(trainpars,'psbasisGrad');
        end
        
        if (isfield(trainpars,'cpbasisGrad'))
            trainpars = rmfield(trainpars,'cpbasisGrad');
        end
        
        if (isfield(trainpars,'kx'))
            trainpars = rmfield(trainpars,'kx');
        end
        
        if (~exist('tag','var') || isempty(tag))
            matfilename = strcat(save_dir,cell_prefix,'_cpl_',form_cell_prefix(cells_touse(N+1:Neff)),'_',file_suffix,'.mat');
        else
            matfilename = strcat(save_dir,tag);
        end
        
        if(isfield(trainpars,'lgrad'))
            trainpars = rmfield(trainpars,'lgrad');
        end
        
        basepars.Nneurons = N;
        basepars.crop_idx = master_crop_idces;
        
        1;
        
        save(matfilename,'filtermode','master_p0','master_pstar','basepars',...
            'stimpars','trainpars','tstim','cell_ids','stim_type', 'fstar', 'gstar','master_crop_idces','cells_totrain','-v7.3');
        
end