function  Simulations = simulate_model(SimPars, CIFcomponents)  
% Checked 2012-12-07   AK Heitman
% Self Contained.. no calls down

%INPUTS: 
%-SimPars (Simulation Parameters)
%   .bin_dt     time in secs of the spike bin
%   .trials     simulations_pertrial
%   .seconds    the total number of seconds in the trial
%   .Coupling   Logical.   is there or isn't there coupling
%   .PostSpike  Logical.   existence of post spike filter
%   .simulations_pertrial    self explanatory
%-CIFcomponents (Conditional Intensity Function)
%   .lcif_mu    basically just mu constant for each spikebin
%   .lcif_kx    filtered stimulus, changes with frames
%               piecewise constant wrt 
%   .lcif_cp    Cell with entries equal to number of raster trials
%               {trial}.couplingterm   = the actual coupling term
%               {trial}.blocknum
%OUTPUTS: Simulations{i_rastertrial, i_simnum}
% .binnedlogical  .. boolean vector of spikes
% .finalcif       .. conditonal intensity including post spike changes
% .poisson        .. boolean indicating whether sim was poisson
%                 .. true for no post spike, false with ps filter
% .blocknum       .. corresponding block of raster and gathered data
% .simnum         .. simulation number for the block
% .sptimes_secs   .. spike times in seconds
% .sprate         .. simply in hz.  the spike rate





% Unpack SimPars
secs      = SimPars.seconds;
sims      = SimPars.simulations_pertrial;
dt        = SimPars.bin_dt;    

% Load trial independent component
lcif_mukx = CIFcomponents.lcif_mu  + CIFcomponents.lcif_kx;

if isfield(CIFcomponents, 'lcif_rec')
    lcif_mukx = lcif_mukx + CIFcomponents.lcif_rec;
end

bins      = length(lcif_mukx);

%%% With Post Spike this is no longer Poisson process
if SimPars.PostSpike
    cif_psgain = CIFcomponents.cif_psgain;
    ps_bins     = length(cif_psgain);
    pctdonepre = 0;
    for i_trial = 1 : SimPars.trials
        if SimPars.Coupling
            lcif0 = lcif_mukx + CIFcomponents.lcif_cp{i_trial}.couplingterm;
        else
            lcif0 = lcif_mukx;
        end
        cif0         = exp(lcif0);         
        for i_Sim = 1 : sims
            cif_ps       = cif0;
            binary_simulation = false(1,bins);
            for i = 1 : bins- ps_bins;
                roll = rand ( 1) ;  
                if roll >  exp(-dt*cif_ps(i));
                    cif_ps(i+1: i + ps_bins) =  cif_ps(i+1: i + ps_bins) .* (cif_psgain);
                    binary_simulation(i)= true;
                end
            end
            positivebins = find(binary_simulation);
            spikesecs = SimPars.seconds * (positivebins / bins);

            Simulations{i_trial, i_Sim}.binnedlogical   = binary_simulation;
            Simulations{i_trial, i_Sim}.finalcif        = cif_ps;
            Simulations{i_trial, i_Sim}.poisson         = false;
            if SimPars.Coupling
                Simulations{i_trial, i_Sim}.raster_trial    = i_trial;  %%% these are 2 related measures         
                Simulations{i_trial, i_Sim}.blocknum        = CIFcomponents.lcif_cp{i_trial}.blocknum;
            end
            Simulations{i_trial, i_Sim}.simnum          = i_Sim;
            Simulations{i_trial, i_Sim}.sptimes_secs    = spikesecs;
            Simulations{i_trial, i_Sim}.sprate          = max(size(spikesecs)) / secs ; 
        end
        pctdone =100* i_trial / SimPars.trials;
        if round(pctdone/20) > round(pctdonepre/20);
            display(sprintf('PercentSimulationsComplete_%d',round(pctdone)));
        end
        pctdonepre = pctdone;          
    end
end


%%% Poisson Version without the Post-Spike Filter
if ~SimPars.PostSpike
    pctdonepre = 0;
    for i_trial = 1 : SimPars.trials
        if SimPars.Coupling
            lcif0 = lcif_mukx + CIFcomponents.lcif_cp{i_trial}.couplingterm;
        else
            lcif0 = lcif_mukx;
        end
        cif0      = exp(lcif0);
        cif       = cif0; 
        for i_Sim = 1 : sims
            binary_simulation = false(1,bins);
            for i = 1 : bins;
                roll = rand ( 1) ;  
                if roll >  exp(-dt*cif(i));
                    binary_simulation(i)= true;
                end
            end
        %    spikes_cp = size(find(binarycp_simulation));
            positivebins = find(binary_simulation);
            spikesecs = SimPars.seconds * (positivebins / bins);

            Simulations{i_trial, i_Sim}.binnedlogical   = binary_simulation;
            Simulations{i_trial, i_Sim}.cif             = cif;
            Simulations{i_trial, i_Sim}.poisson         = true;
            if SimPars.Coupling
                Simulations{i_trial, i_Sim}.raster_trial    = i_trial;  %%% these are 2 related measures         
                Simulations{i_trial, i_Sim}.blocknum        = CIFcomponents.lcif_cp{i_trial}.blocknum;
            end
            Simulations{i_trial, i_Sim}.simnum          = i_Sim;
            Simulations{i_trial, i_Sim}.sptimes_secs    = spikesecs;
            Simulations{i_trial, i_Sim}.sprate          = max(size(spikesecs)) / secs ; 
        end
        pctdone =100* i_trial / SimPars.trials;
        if round(pctdone/20) > round(pctdonepre/20);
            display(sprintf('PercentSimulationsComplete_%d',round(pctdone)));
        end
        pctdonepre = pctdone;        
    end
end
    
            



