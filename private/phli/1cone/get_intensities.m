function intensities = get_intensities(stimstructs, varargin)
% GET_INTENSITIES
% usage: intensities = get_intensities(stimstructs, varargin)
%
% 2011 phli
%

opts = inputParser;
opts.addParamValue('intensity_predicate', []);
opts.addParamValue('with_blanks', true);
opts.parse(varargin{:});
opts = opts.Results;

intensities = [];
for i = 1:length(stimstructs)
    if iscell(stimstructs) stimstruct = stimstructs{i};
    else stimstruct = stimstructs(i); end    
    intensities = [intensities; unique(stimstruct.rgbs(:,1))];
end
intensities = unique(intensities);


if ~opts.with_blanks
    intensities = intensities(intensities ~= 0);
end


if ~isempty(opts.intensity_predicate)
    intensities = select(num2cell(intensities), opts.intensity_predicate);
    intensities = cell2mat(intensities);
end

