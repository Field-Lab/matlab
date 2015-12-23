function [ud_urgbs ud_colors] = ud_rasterurgbs(stimulus, varargin)

opts = inputParser();
opts.addParamValue('colors', []);
opts.addParamValue('cmf', @jet);
opts.parse(varargin{:});
opts = opts.Results;

ncones = stimulus.numcones;
urgb = stimulus.urgb;
uds = urgb.uds;
intensities = urgb.intensities;
singles = urgb.singles;

% Assume this is the max for all singles and U/Ds
maxintensity = max(urgb.absolute_intensities(:));

% Fill in U/Ds
ud_urgbs = cell(ncones);
for i = 1:ncones
    for d = 1:ncones
        if i == d, continue; end
        bool = uds & intensities(i,:) == maxintensity & intensities(d,:) == -maxintensity;
        ud_urgbs{i,d} = find(bool);
    end
end


% Col incr singles
incr_urgbs = cell(ncones,1);
for i = 1:ncones
    incr_urgbs{i} = find(singles & intensities(i,:) == maxintensity);    
end

% Row decr singles
decr_urgbs = cell(1,ncones);
for d = 1:ncones
    decr_urgbs{d} = find(singles & intensities(d,:) == -maxintensity);        
end


% New scheme is to plop the incrs on top, decrs on left
% (reverse for OFF cells; just transpose after)
ud_urgbs = [decr_urgbs; ud_urgbs];
incr_urgbs = [{[]}; incr_urgbs];
ud_urgbs = [incr_urgbs ud_urgbs];


% Create colors matrix?
if nargout > 1
    % Black for U/Ds (don't worry that diagonals are marked)
    ud_colors = repmat({{[0 0 0]}}, ncones+1);
    
    % Colors for rows and columns
    if isempty(opts.colors), opts.colors = opts.cmf(ncones); end
    for i = 1:ncones
        ud_colors{i+1,1} = {opts.colors(i,:)};
        ud_colors{1,i+1} = {opts.colors(i,:)};
    end
end