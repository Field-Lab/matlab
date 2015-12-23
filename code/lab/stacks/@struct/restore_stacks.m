function datarun = restore_stacks(datarun, varargin)
% RESTORE_STACKS    Load saved stack data and merge it into datarun
% usage: datarun = restore_stacks(datarun, opts)
%
% opts: override false  Should loaded data override any existing stacks in
%                       datarun, or should they be combined more
%                       selectively; see MERGE_STACKS for details
%
% Unmatched options are passed along to retrieve_stacks.
%
% See also: MERGE_STACKS, RETRIEVE_STACKS
%
% 2010-10 phli
%

opts = inputParser;
opts.addParamValue('override', false);
opts.KeepUnmatched = true;
opts.parse(varargin{:});
unmatched = opts.Unmatched;
opts = opts.Results;


loaded = retrieve_stacks(datarun, unmatched);

% If we are overriding datarun.stacks, then we're done
if opts.override || ~isfield(datarun, 'stacks')
    datarun.stacks = loaded.stacks;
    return
end

% Otherwise we have a big job: first decide which stacks might be
% inconsistent between datarun and loaded.  Then carefully merge loaded
% into datarun, skipping inconsistencies.
datarun = merge_stacks(datarun, loaded);