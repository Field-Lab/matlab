function weights = calc_stim_rf_weights(datarun, cellid, stimmaps, varargin)
% CALC_STIM_RF_WEIGHTS
% usage: calc_stim_rf_weights(datarun, cellid, stimmaps, [mapindices])
%
% 2012-02 phli
%

opts = inputParser();
opts.addParamValue('rfopts', {});
opts.addParamValue('mapindices', []);
opts.parse(varargin{:});
opts = opts.Results;

rfs{1,1} = get_rf(datarun, cellid, opts.rfopts{:});
rfsize = size(rfs{1,1});

% Kludge to handle single map input nicer
cellout = true;
if ~iscell(stimmaps)
    stimmaps = {stimmaps};
    cellout = false;
end


for i = 1:length(stimmaps)
    stimmap = stimmaps{i};
    scaleup = size(stimmap) ./ rfsize;
    
    % Check whether we already calculated a scaled up version.  Generate it if necessary.
    if any(size(rfs) < scaleup) || isempty(rfs{scaleup(1), scaleup(2)})
        rfs{scaleup(1), scaleup(2)} = matrix_scaled_up(rfs{1,1}, [], struct('scale_x', scaleup(1), 'scale_y', scaleup(2)));
    end
    rf = rfs{scaleup(1), scaleup(2)};

    if isempty(opts.mapindices), opts.mapindices = setdiff(unique(stimmap), 0); end
    weights{i} = zeros(length(opts.mapindices),1);
    for mapindex = opts.mapindices(:)'
        masked = rf(stimmap == mapindex);
        weights{i}(mapindex) = sum(masked(:));
    end
end


% Kludge; see above
if ~cellout, weights = weights{1}; end