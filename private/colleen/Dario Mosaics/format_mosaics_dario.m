% format receptive fields for Dario

% datarun=load_data('2009-12-03-0/data003-nwpca/daat003/daat003.neurons');
datarun=load_data('2008-04-30-1/data021/data021.neurons');
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
    sample = datarun.stas.fits{cell_id};
    xy(i,: ) = sample.mean;
     ang(i,: ) = sample.angle;
     sd(i,: ) = sample.sd;
end
ON.xy = xy;
ON.ang = ang;
ON.sd = sd;
% plot_rf_summaries(datarun, on_cell_ids, ON)
for i = 1: size(ON.sd,1)
    if ON.sd(i,1) < ON.sd(i,2)
        ON.sd(i,:) = [ON.sd(i,2), ON.sd(i,1)];
        if ON.ang(i) < pi/2
            ON.ang(i)= ON.ang(i) + 3*pi/2;
        else
            ON.ang(i) = ON.ang(i)-pi/2;
        end
        if ON.ang(i)>= pi
            ON.ang(i) = ON.ang(i) - pi; 
        end
        
    else
        ON.sd(i,:) = [ON.sd(i,1), ON.sd(i,2)];
         if ON.ang(i)>= pi
            ON.ang(i) = ON.ang(i) - pi; 
        end
    end
end
% plot_rf_summaries(datarun, on_cell_ids, ON)
off_cell_ids = datarun.cell_types{1,2}.cell_ids;
[cell_numbers] = get_cell_indices(datarun, off_cell_ids);
xy = zeros(size(off_cell_ids,2),2);
ang = zeros(size(off_cell_ids,2),1);
sd = zeros(size(off_cell_ids,2),2);

for i = 1:size(off_cell_ids,2)
    cell_id = cell_numbers(i);
    sample = datarun.stas.fits{cell_id};
    xy(i,: ) = sample.mean;
     ang(i,: ) = sample.angle;
     sd(i,: ) = sample.sd;
end
OFF.xy = xy;
OFF.ang = ang;
OFF.sd = sd;

for i = 1: size(OFF.sd,1)
    if OFF.sd(i,1) < OFF.sd(i,2)
        OFF.sd(i,:) = [OFF.sd(i,2), OFF.sd(i,1)];
        if OFF.ang(i) < pi/2
            OFF.ang(i)= OFF.ang(i) + 3*pi/2;
        else
            OFF.ang(i) = OFF.ang(i)-pi/2;
        end
         if OFF.ang(i)>= pi
            OFF.ang(i) = OFF.ang(i) - pi; 
        end
    else
        OFF.sd(i,:) = [OFF.sd(i,1), OFF.sd(i,2)];
        OFF.ang(i) = OFF.ang(i);
         if OFF.ang(i)>= pi
            OFF.ang(i) = OFF.ang(i) - pi; 
        end
    end
end

save('p200912030_d003_RFc', 'ON', 'OFF');