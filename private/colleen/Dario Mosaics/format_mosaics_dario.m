% format receptive fields for Dario

datarun=load_data('2010-10-18-2/data000-nwpca/data000-nwpca.neurons');

% Loading other information
%   Other information you might want including STAs, params, ei, neurons
%   (includes spike times) and polarities
datarun=load_sta(datarun);
datarun=load_params(datarun);
datarun=load_neurons(datarun);
datarun=set_polarities(datarun);

on_cell_ids = datarun.cell_types{1,1}.cell_ids;
[cell_numbers] = get_cell_indices(datarun, on_cell_ids);

xy = zeros(size(on_cell_ids,2),2);
ang = zeros(size(on_cell_ids,2),1);
sd = zeros(size(on_cell_ids,2),2);

for i = 1:size(on_cell_ids,2)
    cell_id = cell_numbers(i);
    sample = datarun.vision.sta_fits{cell_id};
    xy(i,: ) = sample.mean;
     ang(i,: ) = sample.angle;
     sd(i,: ) = sample.sd;
end
ON.xy = xy;
ON.ang = ang;
ON.sd = sd;



off_cell_ids = datarun.cell_types{1,2}.cell_ids;
[cell_numbers] = get_cell_indices(datarun, off_cell_ids);
xy = zeros(size(off_cell_ids,2),2);
ang = zeros(size(off_cell_ids,2),1);
sd = zeros(size(off_cell_ids,2),2);

for i = 1:size(off_cell_ids,2)
    cell_id = cell_numbers(i);
    sample = datarun.vision.sta_fits{cell_id};
    xy(i,: ) = sample.mean;
     ang(i,: ) = sample.angle;
     sd(i,: ) = sample.sd;
end
OFF.xy = xy;
OFF.ang = ang;
OFF.sd = sd;
save('p201010182_d000_RFc', 'ON', 'OFF');