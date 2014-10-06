function datarun = load_txt_cell_types(datarun, txtpath, varargin)
% STRUCT/LOAD_TXT_CELL_TYPES   Datarun wrapper for LOAD_TXT_CELL_TYPES
% usage: datarun = load_txt_cell_types(datarun, txtpath, opts)
%
% See also: LOAD_TXT_CELL_TYPES
%

opts = inputParser();
opts.addParamValue('verbose', true);
opts.KeepUnmatched = true;
opts.parse(varargin{:});
unmatched = opts.Unmatched;
opts = opts.Results;

% Try to find TXTPATH in different ways under RRS_PREFIX
if ~exist(txtpath, 'file') && ~strcmp(txtpath(1),'/') && ~is_windows_root_path(txtpath)
    path1 = sprintf('%s.%s',     datarun.names.rrs_prefix, txtpath);
    path2 = sprintf('%s%s',      datarun.names.rrs_prefix, txtpath);
    path3 = sprintf('%s.%s.txt', datarun.names.rrs_prefix, txtpath);
    path4 = sprintf('%s%s.txt',  datarun.names.rrs_prefix, txtpath);
    
    for testpath = {path1 path2 path3 path4}
        if exist(testpath{1}, 'file')
            txtpath = testpath{1};
            break
        end
    end
end

if ~exist(txtpath, 'file')
    warning('STRUCT/LOAD_TXT_CELL_TYPES: could not find classification file %s', txtpath);
    return
end


% Have to turn unmatched back into a cell array to allow sending to the
% non-STRUCT overloaded version of LOAD_TXT_CELL_TYPES...
load_celltype_args = namedstruct2cell(unmatched);

celltypes = load_txt_cell_types(txtpath, load_celltype_args{:});
datarun = load_celltypes_into_datarun(datarun, celltypes, opts.verbose);