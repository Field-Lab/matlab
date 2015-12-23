function datarun = calc_csd(datarun, varargin)
% CALC_CSD      Calculate CSDs based on EIs and array info and save into datarun
%
% It is important not to have floating point inaccuracy when calculating
% the Delaunay Triangulation of the electrode positions.  Therefore, the
% electrode positions are rounded off, which assumes that they are
% generally whole numbers.  Will try to check that this is not violated and
% throw a warning if it is.
%
% The sparse CSD multiplication matrix returned by CSDMATRIX has zeros in
% positions for bad electrodes and electrodes along the edge for which the
% divergence cannot be properly calculated.  If the zeros will be
% confusing, should add option to replace them with NaNs.  But for general
% EI plotting having zeros in these positions is probably okay.
%
% 2013-05, phli
%

opts = inputParser();
opts.addParamValue('positions', []);
opts.addParamValue('badelectrodes', NaN);
opts.addParamValue('cellspec', 'all');
opts.parse(varargin{:});
opts = opts.Results;


csdcellnums = get_cell_indices(datarun, opts.cellspec);


% Load EIs as needed
if ~isfield(datarun, 'ei') || ~isfield(datarun.ei, 'eis') || isempty(datarun.ei.eis)
    missingspec = 'all';
else
    missingnums = cellfun(@isempty, datarun.ei.eis);
    loadnums = intersect(csdcellnums, missingnums);
    missingspec = datarun.cell_ids(loadnums);
end
datarun = load_ei(datarun, missingspec);
if isempty(opts.positions), opts.positions = datarun.ei.position; end


% Check that electrode positions can be safely rounded off
positions = round(opts.positions);
err = max(max(abs(positions - opts.positions)));
if err > 1e-6
    warning('CALC_CSD:ROUNDING', 'Electrode positions should be upscaled to allow rounding off');
end


% Load disconnected if nothing was given
if isnan(opts.badelectrodes)
    opts.badelectrodes = [];
    if isfield(datarun.ei, 'disconnected')
        opts.badelectrodes = datarun.ei.disconnected;
    end
end


% Get CSD matrix (the unclipped Delaunay Triagulation could be cached)
csdM = csdmatrix(round(datarun.ei.position), 'badelectrodes', opts.badelectrodes);


% Calculate CSDs
ncsds = length(csdcellnums);
datarun.ei.csds = cell(ncsds,1);
for i = 1:ncsds
    cellnum = csdcellnums(i);
    datarun.ei.csds{cellnum} = csdM * datarun.ei.eis{cellnum};
end