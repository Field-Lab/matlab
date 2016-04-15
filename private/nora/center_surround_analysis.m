function center_surround_analysis(varargin)

p = inputParser;
p.addParameter('sigmas', 2)
p.addParameter('center', 'k')
p.addParameter('surround', 'k')
p.addParameter('full', 'k')
p.parse(varargin{:});

n_sigmas = length(p.Results.sigmas);
n_cells = size(p.Results.center{1},2);
for i_sigma = 1:n_sigmas
    for i_cell = 1:n_cells
        if iscell(center) && iscell(full)
            PSTH_full = IDP_plot_PSTH(full, i_cell, 0, 1/120, 10);
            PSTH_center = IDP_plot_PSTH(center, i_cell, 0, 1/120, 10);
            NSEM_Corr(i_cell,i_sigma) = err(PSTH_full, PSTH_center);
        end
        if iscell(surround)
            PSTH_surround = IDP_plot_PSTH(surround, i_cell, 0, 1/120, 10);
            surround_struct(i_cell,i_sigma) = std(PSTH_surround);
        end
    end
end

end