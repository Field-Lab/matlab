% Function to simulate the GLM using paramters specified by (b,K,PS,CP)
% basepars) driven by a stimulus X for ntrials

% Returns the binary matrix of spikes and the log of the resulting CIFs (time x neuron x trial)

function [spikes lcifs hist_term kx STAwindows] = simulate_model_nopars(basepars,b,K,PS,CP,X,Deff,dt,ntrials,opt_pars)


1;

% Preprocess stimulus if needed
if (isfield(basepars,'frontend'))
    X = basepars.frontend(X);
end

if (size(PS,1) < size(PS,2))
    PS = PS';
    fprintf('Transposing postspike filter vector: dimensions are now %d by %d\n',size(PS,1),size(PS,2));
end

% Number of time bins in spike train/cif
1;
spT = basepars.maxt * basepars.fac;

% Initalize output
lcifs = zeros(spT,basepars.Nneurons,ntrials);
hist_term = zeros(spT,basepars.Nneurons,ntrials);
spikes = (logical(false(spT,basepars.Nneurons,ntrials)));

% Setup a stimpars struct just to filter it
stimpars.x = X;
stimpars.dt = dt*basepars.fac;

% Filter the stimulus
if (~exist('opt_pars','var'))
    kx = filter_stimulus_direct(K,X,basepars,stimpars.dt);
    %kx = filterstimulus_train2_nonsep(p,basepars,stimpars,[]);
    1;
elseif (isfield(opt_pars,'X'))
1;
    kx = filter_stimulus_direct(K,opt_pars.X,basepars,stimpars.dt);
    1;
    fprintf('Using supplied stimulus!\n');
end

N = basepars.Nneurons;
Neff = N + size(Deff,2);

exflag = 1;
1;
MAX_ITER = 3;

% Truncate the postspike filter when possible
pstol = 10^(-10);

lastnonzero = find(abs(flipud(PS)) > pstol,1);
if (isempty(lastnonzero))
    fprintf('PS filter seems to be zero...');
    cutoff = 1;
else
    cutoff = size(PS,1) - lastnonzero;
end
if (~isempty(cutoff))
    basepars.Mhist = cutoff;
    PS = PS(1:cutoff,:);
    fprintf('Truncated postspike filter at %d samples.\n',cutoff);
end
    


sta_flag = exist('opt_pars','var') && isfield(opt_pars,'sta_windows') && nargout > 3;


if (sta_flag)
        STAwindows = zeros(basepars.n,basepars.Mk,opt_pars.nstawindows);
        NSP = zeros(opt_pars.nstawindows,1);
else
      STAwindows = [];
end

1;

disppercent(-inf,'Simulating trials.');
for j=1:ntrials

    % set a seed for exprnd
    %rand('twister',ntrials*j);
    
    counter = 1;

    if (mod(j,10) == 0)
        disppercent(j/ntrials);    
        %fprintf('Completed %d trials.\n',j);
    end
    
    % mystream =  RandStream.create('mt19937ar','Seed',100*clock);
    
    while(exflag && counter < MAX_ITER)
        
        if (counter > 1)
            1;
            %mystream = RandStream.create('mt19937ar','Seed',MAX_ITER*j+counter);
            fprintf('Restarting trial %d: iteration %d\n',j,counter);
        end
        
        
        if (exist('opt_pars','var') && ~isfield(opt_pars,'X'))
   
           %Regenerate random stimulus
           switch (opt_pars.mode)
              
               case 'gwn'
                   stimpars.x = opt_pars.mu + opt_pars.sig.*randn(basepars.n,basepars.maxt);
               
               case 'switch'
                   
                   eff_end = opt_pars.period * floor(basepars.maxt/opt_pars.period);
                   numperiods = eff_end/opt_pars.period;
                   stimpars.x = zeros(basepars.n,basepars.maxt);
                   for k=1:numperiods
                       idx = (k-1)*opt_pars.period+1:k*opt_pars.period;
                   
                       switch(opt_pars.noisetype)
                           case 'gwn'
                               stimpars.x(:,idx) = [(opt_pars.mu1+opt_pars.sig1.*randn(basepars.n,floor(opt_pars.period/2))) (opt_pars.mu2+opt_pars.sig2.*randn(basepars.n,ceil(opt_pars.period/2)))];
                           case 'bwn'
                               1;
                               stimpars.x(:,idx) = [(opt_pars.mu1-opt_pars.sig1 + 2.*opt_pars.sig1.*(rand(basepars.n,floor(opt_pars.period/2))>0.5)) (opt_pars.mu2-opt_pars.sig2 + 2.*opt_pars.sig2.*(rand(basepars.n,ceil(opt_pars.period/2))>0.5))];
                       end
                   end
                   1;
                   % Fill in the remainder with 1st config
                   switch(opt_pars.noisetype)
                       case 'gwn'
                        stimpars.x(:,eff_end+1:end) = opt_pars.mu1 + opt_pars.sig1.*randn(basepars.n,basepars.maxt-eff_end);
                       case 'bwn'
                        stimpars.x(:,eff_end+1:end) = opt_pars.mu1-opt_pars.sig1 + 2.*opt_pars.sig1.*(rand(basepars.n,basepars.maxt-eff_end)>0.5);
                   end
           end
           kx = filter_stimulus_direct(K,stimpars.x,basepars,stimpars.dt);
          
           %kx = filterstimulus_train2_nonsep(p,basepars,stimpars,[]);
            1;
        end
        
        1;
        if (isempty(Deff))
             [exflag sp lcifs(:,:,j) hist_term(:,:,j)] = simulate_trains_nopars(basepars,dt,kx,b,[],PS,CP);
        else
            [exflag sp lcifs(:,:,j) hist_term(:,:,j)] = simulate_trains_nopars(basepars,dt,kx,b,Deff(:,:,j),PS,CP);
        end
        counter = counter+1;
    end
    
   
    if (exflag)
        spikes = [];
        fprintf('Failed to simulate with this model.\n');
        spikes = [];
        lcifs = [];
        hist_term = [];
        kx = [];
        return;
    end        
    exflag = 1;
    
    
    
    [I J] = find(sp);
    for k=1:length(I)
        spikes(I(k),J(k),j) = true;
    end
    
    1;
    
    if (sta_flag)
        % Compute STA using specified time windows
        
        if (sum(NSP < opt_pars.max_spikes) == 0)
            fprintf('All windows have been computed!Returning out.\n')
            return;
        end
        
        
        for k=1:opt_pars.nstawindows
            
            if (NSP(k) < opt_pars.max_spikes)
                
                for r=2:numperiods
                    frame_offset = (r-1)*opt_pars.period + opt_pars.sta_windows(k);
                    sp_window = dt.*I(I> basepars.fac*frame_offset & I< basepars.fac*(frame_offset + opt_pars.windowsize)) - dt.*(basepars.fac*frame_offset);
                    %fprintf('Computing STA for trial %d frameset=[%d,%d]\n',j,frame_offset+1,frame_offset+opt_pars.windowsize);
                    [fsta nsp] = fast_sta(stimpars.x(:,frame_offset+1:frame_offset+opt_pars.windowsize),sp_window,basepars.Mk,stimpars.dt);
                    if (nsp > 0)
                        STAwindows(:,:,k) = STAwindows(:,:,k) + nsp*fsta;
                    end
                    NSP(k) = NSP(k) + nsp;
                    fprintf('Window %d has accumulated %d spikes.\n',k,NSP(k));
                end
                
            else
                fprintf('Completed STA for window %d. Skipping...\n',k);
            end
            
        end
    end
    
    
end
elapsedTime = disppercent(inf);

if (sta_flag)
    for k=1:opt_pars.nstawindows
        STAwindows(:,:,k) = STAwindows(:,:,k)./NSP(k);
    end
end

