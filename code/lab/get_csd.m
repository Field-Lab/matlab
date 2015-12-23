function csd = get_csd(datarun, cellid)
cellnum = get_cell_indices(datarun, cellid);

if isfield(datarun.ei, 'csds') && length(datarun.ei.csds) >= cellnum && ~isempty(datarun.ei.csds{cellnum})
    csd = datarun.ei.csds{cellnum};
    return;
end

datarun = calc_csd(datarun, 'cellspec', cellid);
csd = datarun.ei.csds{cellnum};