function datarun = stimulus_from_globals(datarun)

if ~isfield(datarun, 'globals') || isempty(datarun.globals)
    datarun = load_globals(datarun);
end

datarun.stimulus = stimulus_from_globals(datarun.globals);