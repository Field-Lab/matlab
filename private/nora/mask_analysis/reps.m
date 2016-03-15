%% Dataruns
regA = 'data005';
center = 'data007';
LES = 'data008';
surround = 'data009';
surrLES = 'data010';
regB = 'data011';
cells = 'Maskin';
classification = 'data002';

%% classification run 
class_datarun = load_data(['/Volumes/Analysis/2015-12-18-2/streamed/' classification '/' classification]);
class_datarun = load_params(class_datarun);
class_datarun = load_neurons(class_datarun);
% check stability
datarun = load_data(['/Volumes/Analysis/2015-12-18-2/map-from-data002-stream/' regA '/' regA]);
datarun = load_neurons(datarun);
regA_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells, 'datarun_class', class_datarun);
datarun = load_data(['/Volumes/Analysis/2015-12-18-2/map-from-data002-stream/' regB '/' regB]);
datarun = load_neurons(datarun);
regB_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells, 'datarun_class', class_datarun);

%% full A vs full B
for i_cell = 1:6
   regA_spikes = floor(120*cell2mat(regA_data.testspikes(:,i_cell))); 
   regB_spikes = floor(120*cell2mat(regB_data.testspikes(:,i_cell))); 
   for i = 1:1200
       regA_rate(i) = sum(regA_spikes == i);
       regB_rate(i) = sum(regB_spikes == i);
   end
   plot(regA_rate)
   hold on; plot(regB_rate); hold off
   pause()
end

%%
datarun = load_data(['/Volumes/Analysis/2015-12-18-2/map-from-data002-stream/' center '/' center]);
datarun = load_neurons(datarun);
center_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells, 'datarun_class', class_datarun);
datarun = load_data(['/Volumes/Analysis/2015-12-18-2/map-from-data002-stream/' surround '/' surround]);
datarun = load_neurons(datarun);
surround_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells, 'datarun_class', class_datarun);

%% center vs surround
for i_cell = 1:6
   center_spikes = floor(120*cell2mat(center_data.testspikes(:,i_cell))); 
   surround_spikes = floor(120*cell2mat(surround_data.testspikes(:,i_cell))); 
   for i = 1:1200
       center_rate(i) = sum(center_spikes == i);
       surround_rate(i) = sum(surround_spikes == i);
   end
   plot(center_rate)
   hold on; plot(surround_rate); hold off;
   pause()
end

%% center vs full
for i_cell = 1:6
    regB_spikes = floor(120*cell2mat(regB_data.testspikes(:,i_cell)));
    for i = 1:1200
        regB_rate(i) = sum(regB_spikes == i);
    end
    center_spikes = floor(120*cell2mat(center_data.testspikes(:,i_cell)));
    surround_spikes = floor(120*cell2mat(surround_data.testspikes(:,i_cell)));
    for i = 1:1200
        center_rate(i) = sum(center_spikes == i);
        surround_rate(i) = sum(surround_spikes == i);
    end
    plot(center_rate)
    hold on; plot(regB_rate); hold off;
    title(i_cell)
    pause()
end  

%%
datarun = load_data(['/Volumes/Analysis/2015-12-18-2/map-from-data002-stream/' surrLES '/' surrLES]);
datarun = load_neurons(datarun);
surrLES_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells, 'datarun_class', class_datarun);

%% center vs full
for i_cell = 1:6
    regB_spikes = floor(120*cell2mat(regB_data.testspikes(:,i_cell)));
    for i = 1:1200
        regB_rate(i) = sum(regB_spikes == i);
    end
    surLES_spikes = floor(120*cell2mat(surrLES_data.testspikes(:,i_cell)));
    for i = 1:1200
        surrLES_rate(i) = sum(surLES_spikes == i);
    end
    plot(regB_rate)
    hold on; plot(surrLES_rate); hold off;
    title(i_cell)
    pause()
end  

