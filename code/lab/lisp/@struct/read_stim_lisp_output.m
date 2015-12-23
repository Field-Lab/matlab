function datarun = read_stim_lisp_output(datarun, mappath, force, parsepoly)
% READ_STIM_LISP_OUTPUT     Wrap for datarun.
%
% 2012-06 phli
%

if nargin < 2, mappath = []; end
if nargin < 3, force = false; end
if nargin < 4, parsepoly = true; end

if isfield(datarun, 'stimulus') && ~isempty(datarun.stimulus) && ~force
    return
end


spec = parse_rrs_prefix(datarun);
datarun.stimulus = read_stim_lisp_output(spec.piece_fullname, sprintf('s%.2d', spec.run_num), mappath, parsepoly);