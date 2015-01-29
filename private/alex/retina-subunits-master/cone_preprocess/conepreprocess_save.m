function datarun = conepreprocess_save(inrun, varargin)

opts = inputParser();
opts.addParamValue('date',[]);
opts.addParamValue('cone_data_ind', []);
opts.addParamValue('cell_types', {1,2,3,4,5});
opts.parse(varargin{:});
opts = opts.Results;

% Just in case it's not already here
if ~isfield(inrun, 'cell_types')
    inrun = load_params(inrun);
end


% Jeremy wants RFs for sure
if all(cellfun(@isempty, inrun.stas.rfs))
    inrun = get_sta_summaries(inrun, opts.cell_types, 'keep_stas', false, 'marks_params', struct('robust_std_method', 6));
end

% Reduce filesize by cutting out some fields or converting to smaller datatypes
datarun = conepreprocess_strip(inrun);


% Save it
parsed = parse_rrs_prefix(inrun);
savepath = [conepreprocess_path parsed.piece_fullname];
% mkdir(savepath);
savepath =  [savepath '/subunits/' parsed.run_full_name];
if ~isdir(savepath)
    mkdir(savepath);
end
if ischar(opts.cone_data_ind)
    savepath = [savepath '/' opts.cone_data_ind '/' opts.date '/' parsed.run_full_name];
    mkdir(savepath);
end
savepath =  [savepath '/conepreprocess.mat'];

if exist(savepath, 'file')
    overwrite = questdlg(['File ' savepath ' exists; overwrite?'], 'Overwrite existing preprocessing?', 'Okay', 'Cancel', 'Cancel');

    % If aborted saving, check whether we're going to lose all the work and
    % maybe dump to be safe.
    if strcmp(overwrite, 'Cancel')
        if nargout == 0
            warning('Dumping datarun to variable CONEPREPROCESS_DUMP');
            assignin('base', 'CONEPREPROCESS_DUMP', datarun);
        end
        return
    end
    
end

save(savepath, 'datarun');


if nargout == 0
    clear datarun;
end
