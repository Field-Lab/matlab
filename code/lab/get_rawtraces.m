function traces = get_rawtraces(data_spec, centers, nrpoints, varargin)
% GET_RAWTRACES    Pull windowed traces from raw data
% usage: traces = get_rawtraces(data_spec, centers, nrpoints, [opts])
%
% inputs: data_spec    See get_rawdatafile for details
%         centers      The center points for the data windows
%         nrpoints     The number of points right of center to include
%
% options: nlpoints    0   The number of points left of center to include
%          electrodes  []  Which electrodes to include in output.  If left
%                          blank, defaults to all except for channel 0
%                          which is typically the trigger channel.
%
% output: Traces is an IxJxK matrix of traces as column vectors.  I =
% nrpoints + nleftpoints + 1.  J is the number of electrodes.  K is the
% number of center points given.
%
% 2010-06 phli
%

opts = inputParser;
opts.addParamValue('nlpoints', 0);
opts.addParamValue('electrodes', []);
opts.parse(varargin{:});
opts = opts.Results;


rdf = get_rawdatafile(data_spec);
if isempty(opts.electrodes)
    % Assume we want all electrodes, dropping trigger channel 0
    opts.electrodes = (2:length(rdf.getData(1,1))) - 1;
end


num_samples = opts.nlpoints + nrpoints + 1;
num_elec = length(opts.electrodes);
traces = zeros(num_samples, num_elec, length(centers));
for i = 1:length(centers)
    start_point = centers(i) - opts.nlpoints;
    traces_block = rdf.getData(start_point, num_samples);
    traces(:,:,i) = traces_block(:, opts.electrodes + 1);
end