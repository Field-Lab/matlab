% Function to simulate a network of GLM neurons for ntrials trials
% Chaitu Ekanadham 07/27/2009

% Arguments:

%  basepars - base parameters - the field basepars.n will NOT be used
%  p - cell array containing the parameters for each neuron in vector format
%  cpl_mat - N x N logical matrix where jth row specifies which neurons affect cell j via coupling filters
%  X - stimulus matrix (space x time) - crop indices index the rows of this matrix
%  Deff - any spikes you want to condition on 
%  dt - spike train resolution at which to simulate - stimulus resolution is assumed to be dt*basepars.fac
%  ntrials - number of trials to simulate
%  opt_pars - special stimulus generation (documented elsewhere).

% Returns the binary matrix of spikes and the log of the resulting CIFs (time x neuron x trial)
% Optionally returns the filtered stimulus (time x neuron x trial)

function [spikes lcifs kx] = simulate_network_multisize(basepars,p,cpl_mat,X,Deff,dt,ntrials,opt_pars)

basepars.maxt = size(X,2);

% Number of time bins in spike train/cif
spT = basepars.maxt * basepars.fac;

% Initalize output
lcifs = zeros(spT,basepars.Nneurons,ntrials);
spikes = (logical(false(spT,basepars.Nneurons,ntrials)));

% Setup a stimpars struct just to filter it
stimpars.x = X;
stimpars.dt = dt*basepars.fac;

% Filter the stimulus
if (~exist('opt_pars','var'))
    switch(basepars.filtermode)
        case 'sep_raw'
            1;
            kx = filterstimulus_train_cellmode(p,basepars,stimpars);
        case 'nonsep'
            kx = filterstimulus_train_nonsep_cellmode(cell2mat(p),basepars,stimpars,[]);
    end
end

N = basepars.Nneurons;
Neff = N + size(Deff,2);


% Recover base rate, postspike filter  and coupling filter coefficients for these neurons
[b ps_weights cp_weights] = recover_multisize_pars(basepars,p,cpl_mat);



exflag = 1;
1;
MAX_ITER = 10;
disppercent(-inf,'Starting simulation.');
for j=1:ntrials
    
    % set a seed for exprnd
    rand('twister',j);
    
    counter = 1;

    %if (mod(j,10) == 0)
        disppercent(j/ntrials);   
        %fprintf('Completed %d trials.\n',j);
    %end
        
    while(exflag && counter < MAX_ITER)
        
        if (counter > 1)
            fprintf('Restarting trial %d: iteration %d\n',j,counter);
        end
        
        
        if (exist('opt_pars','var'))
           %Regenerate random stimulus
           switch (opt_pars.mode)
              
               case 'gwn'
                   sitmpars.x = opt_pars.mu + opt_pars.sig.*randn(basepars.n,basepars.maxt);
               
               case 'switch'
                   stimpars.x = [(opt_pars.mu1+opt_pars.sig1.*randn(basepars.n,basepars.maxt/2)) (opt_pars.mu2+opt_pars.sig2.*randn(basepars.n,basepars.maxt/2))];
                   
           end
           switch(basepars.filtermode)
               case 'sep_raw'
                    kx = filterstimulus_train_cellmode(p,basepars,stimpars,[]);
               case 'nonsep'
                    kx = filterstimulus_train_nonsep_cellmode(p,basepars,stimpars,[]);
           end
            1;
        end
        
        1;
        
        if (isempty(Deff))
             [exflag sp lcifs(:,:,j)] = simulate_trains_fast4(basepars,dt,kx,b,[],ps_weights,cp_weights,cpl_mat);
        else
            [exflag sp lcifs(:,:,j)] = simulate_trains_fast4(basepars,dt,kx,b,Deff(:,:,j),ps_weights,cp_weights,cpl_mat);
        end
        counter = counter+1;
    end
    
    if (exflag)
        spikes = [];
        fprintf('Failed to simulate with this model.\n');
        return;
    end        
    exflag = 1;
    
    
    
    [I J] = find(sp);
    for k=1:length(I)
        spikes(I(k),J(k),j) = true;
    end
end
disppercent(inf);