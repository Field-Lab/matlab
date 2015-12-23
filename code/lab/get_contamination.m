function datarun = get_contamination(datarun,cell_spec,varargin)
% CONTAMINATION   Calculate contamination of spike train
%
% usage:  c = contamination(...)
%
% Calculate contamination of spike train. Optionally return the 
% number of refractory violations.
%
% shlens 2006-07-17
% greschner


% SET UP OPTIONAL ARGUMENTS
    p = inputParser;
    p.addParamValue('field_name', 'contamination');% 
    p.addParamValue('acf_range', [0.5 1.2] * 1e-3);% range of ACF considered contaminated  
    % parse inputs
    p.parse(varargin{:});
    params = p.Results;

    
%get field       
if isfield(datarun, params.field_name);
    conta=datarun.(params.field_name);
else
    conta=zeros(length(length(datarun.cell_ids)),1);
end


cell_indices = get_cell_indices(datarun,cell_spec);

delta = 1e3 * (params.acf_range(2) - params.acf_range(1));
options = struct('dt',0.05e-3,'offset',50e-3,'max_isi',0);

for i=1:length(cell_indices)
    
    % calculate ACF (50 ms offset, 0.05 ms dt
    [acf, time]=compute_ccf(datarun.spikes{cell_indices(i)}, datarun.spikes{cell_indices(i)},options);

    % convert ACF to raw counts
    acf_counts = acf * (length(datarun.spikes{cell_indices(i)}) * options.dt);

    % determine relevant time
    violation_time = find((time <= params.acf_range(2)) & (time >= params.acf_range(1)));

    % count violations
    violations = sum(acf(violation_time));

    % duration in milliseconds
    duration = datarun.spikes{cell_indices(i)}(end) * 1e3;

    % number of spikes
    N = length(datarun.spikes{cell_indices(i)});

    % calculate contamination
    conta(cell_indices(i)) = violations * duration / (delta * N * N);

end


datarun.(params.field_name)=conta; 


