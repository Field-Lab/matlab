%% Clean sparse matrix multiplication method with Delaunay compensation for bad electrode positions
csdM = csdmatrix(ai.positions, 'badelectrodes', ai.disconnected);
csd = csdM * ei;

% Get all of them; very fast once EIs are loaded
datarun = calc_csd(datarun);


%% Visualization
% Normal EI plot
plot_ei_scroll(datarun, datarun.cell_ids(91));

% CSD plot
plot_ei_scroll(datarun, datarun.cell_ids(91), 'modify_ei', @(ei, datarun, cell_id) (ei2csd(ei, neighbor_struct)));

% Plot by maximums
plot_ei_(datarun.ei.eis{91}, datarun.ei.position, 0, 'neg_color', [0 0 1]);
plot_ei_(ei2csd(datarun.ei.eis{91}, neighbor_struct), datarun.ei.position, 0, 'neg_color', [0 0 1]);


%% Deprecated kludgy method

% I didn't come back to this for a while, so not sure if this is really
% even the right way to do it, but it seems to work...

% Get neighbor_struct for calculating CSDs
for i = 1:519
    neighbor_struct(i) = get_lazyhex_ei_neighbors(i, setdiff(get_ei_neighbors(i, 519), i), datarun.ei.position);
    i;
end
