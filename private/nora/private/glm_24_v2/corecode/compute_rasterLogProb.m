% AK HEITMAN
% Reviewed 2012-12-11

% LogProb of each Retina sp train given proposed model (CIFcomponents)

%INPUT: RasterCell , CIFcomponents , duration
%-RasterCell    Raster trials number of entries
%  {i}.spikes   spike times of ret on trial i in seconds
%-CIFcomponents (Conditional Intensity Function)
%   .lcif_mu    basically just mu constant for each spikebin
%   .lcif_kx    filtered stimulus, changes with frames
%               piecewise constant wrt 
%   .lcif_cp    Cell with entries equal to number of raster trials
%               {trial}.couplingterm   = the actual coupling term
%               {trial}.blocknum
%OUTPUT: RASTER_LOGPROB
% .perBin       vector ofthe logprob of each bin, collapse out trials
%               each trial of Bin independent event..simple summation
% .perTrial     cell with trials entries
%               each Bin of each trial ind event.. simple summation
% ._avgProbability
%               additive avg of LogProbs then exponentiated
%               good measure of avg probability of each individiual bin
%Raster only has spike times, all computations based upon CIFcomponents

function Raster_LogProb = compute_rasterLogProb(RasterCell , CIFcomponents, duration)


%% Unpack the Model Components
lcif_mu = CIFcomponents.lcif_mu;  lcif_kx = CIFcomponents.lcif_kx; 
Coupling = false ;     PostSpike = false ; Rectification = false;
if isfield(CIFcomponents, 'lcif_cp')   ,    lcif_cp = CIFcomponents.lcif_cp;     Coupling      = true; end
if isfield(CIFcomponents, 'cif_psgain'), cif_psgain = CIFcomponents.cif_psgain;  PostSpike     = true; end
if isfield(CIFcomponents, 'lcif_rec'),   lcif_rec   = CIFcomponents.lcif_rec;    Rectification = true; end
trials = max(size(RasterCell));
bins   = length(lcif_mu);
dt     = duration / bins;

%%  Solve for the Final Conditional Intensity Function of the Trial
%   Use Retina spike times for PS and CP filter
%   Use lcif_mu and lcif_kx already solved for by the GLM

CIFmatrix         = zeros(bins,trials);
LogicalSpike      = false(bins,trials); 
AntiLogicalSpike  = true (bins,trials);
for i_trial = 1:trials
    if  Coupling, lcif0 = lcif_mu + lcif_kx + lcif_cp{i_trial}.couplingterm; end
    if ~Coupling, lcif0 = lcif_mu + lcif_kx; end
    
    if Rectification, lcif0 = lcif0 + lcif_rec;   end
    
    % Modify the CIF via PS Filter after seeing spike
    SpikeBins = ceil( RasterCell{i_trial}.spikes / dt ); 
    cif0 = exp(lcif0);
    cif  = cif0;
    if PostSpike 
        psbins    = length(cif_psgain);
        for spnum = 1 : length(SpikeBins)
            cifchange = ( SpikeBins(spnum) + 1) : ( SpikeBins(spnum) + psbins) ;
            if max(cifchange)>=bins
                cifchange = (SpikeBins(spnum) + 1) : bins  ;   
            end
            psind = 1:length(cifchange);
            if (SpikeBins(spnum) + 1) >= bins
                cifchange = [];
                psind = [];
            end           
            cif( cifchange ) = cif( cifchange) .*cif_psgain(psind);
        end
    end
    CIFmatrix(:, i_trial)               = cif;
    LogicalSpike(SpikeBins,i_trial)     = true;
    AntiLogicalSpike(SpikeBins,i_trial) = false;
    
end

%% Not Efficient but Mathematically Clean.. Calculation of Log Prob in full

%%%%%%%%%%%%%% COMPUTE LOGPROB OF SEPERATE SPIKE AND NONSP EVENTS
nospike_logprob_MATRIX  = (-dt * CIFmatrix);
spike_logprob_MATRIX    = log( ones(size(CIFmatrix)) - exp(1).^(-dt*CIFmatrix) );
Raster_LogProb_Spikes   = LogicalSpike .* spike_logprob_MATRIX;
Raster_LogProb_NoSpikes = AntiLogicalSpike .* nospike_logprob_MATRIX; 
%%%%%%%%%%%%%% LOGPROB OF EACH APPROPRIATE UNIT 
LogProb_BinTrial        = Raster_LogProb_Spikes +  Raster_LogProb_NoSpikes;
LogProb_perBin          = sum(LogProb_BinTrial,2);
LogProb_pertrial        = sum(LogProb_BinTrial,1);  
%%%%%%%%%%%%%% FILL IN OUTPUT
%%%%%%%%%%%%%% COMPUTE "avgProbability" OF EACH BIN IN UNIT
%%%%%%%%%%%%%% "avgProbability" IS MONOTONIC TRANSFORM OF THE LOG_PROB
LogProb.avgProbability_note = sprintf('Fair Monotonic Transformation...Additive Avg Log Prob run through an Exponential');
LogProb.perBin_note         = sprintf('Look at each trial of given Bin as Independent.  Given CIF of those independent trials, what was likelihood of that string of outcomes');
LogProb.perBin              = LogProb_perBin;
LogProb.perBin_avgProbability = exp(LogProb_perBin/trials);
LogProb.perTrial_note = sprintf('LogProb of each Trial of Raster. Cell Format');
LogProb.perTrial      = cell(1,trials);
for i_trial = 1:trials
    LogProb.perTrial{i_trial}.block          = RasterCell{i_trial}.block;
    LogProb.perTrial{i_trial}.logprob        = LogProb_pertrial(:,i_trial);
    LogProb.perTrial{i_trial}.avgProbability = exp(LogProb_pertrial(:,i_trial) / bins);
end
LogProb.perBinperTrial_note = sprintf('LogProb of each Bin of each Trial, Rows = Bins, Columns = Trials'); 
LogProb.perBinperTrial = LogProb_BinTrial;

Raster_LogProb = LogProb;
end





   % RECALL:   %%--nospikeprob=exp(-dt*cif),spikeprob=(1-exp(-dt*cif))--%%
   % nospike_logprob = (-dt * cif);   
   % spike_logprob   = log(1 - exp(-dt*cif) );
   % AllBins = 1:bins;
   % NoSpikeBins =  setdiff(AllBins , SpikeBins);   
   % LogProb{i_trial}  = sum(spike_logprob(SpikeBins)) +
   % sum(nospike_logprob(NoSpikeBins)); 
