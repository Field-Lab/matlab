class = load_data('/Volumes/Analysis/2015-11-09-8/data029-data033/data029/data029');
NSEM = load_data('/Volumes/Analysis/2015-11-09-8/data029-data033/data032/data032');
EH = load_data('/Volumes/Analysis/2015-11-09-8/data029-data033/data033/data033');
class = load_params(class);
NSEM = load_neurons(NSEM);
EH = load_neurons(EH);
%%
prepped_data_NSEM = interleaved_data_prep(NSEM, 1200, 40, 'cell_spec', 'On Parasol', 'datarun_class', class);
prepped_data_EH = interleaved_data_prep(EH, 1200, 40, 'cell_spec', 'On Parasol', 'datarun_class', class);

%%

bin = 1; % per frame
NSEM_PSTH = zeros(1200*bin,1);
EH_PSTH = zeros(1200*bin,1);

for i_cell = 1:145
NSEM_spikes = floor(cell2mat(prepped_data_NSEM.testspikes(:,i_cell))*120*bin);
EH_spikes = floor(cell2mat(prepped_data_EH.testspikes(:,i_cell))*120*bin);

for i = 1:(120*10*bin)
    NSEM_PSTH(i) = sum(NSEM_spikes ==i);
    EH_PSTH(i) = sum(EH_spikes ==i);
end

plot(conv(NSEM_PSTH, gausswin(10), 'same'))
hold on
plot(conv(EH_PSTH, gausswin(10), 'same'))
hold off
pause()
end

%%
figure;
for i_cell = 1:145
    hold on
    for i = 1:40
        plot(prepped_data_NSEM.testspikes{i, i_cell}, i*ones(length(prepped_data_NSEM.testspikes{i, i_cell})), 'k.')
        plot(prepped_data_EH.testspikes{i, i_cell}, i*ones(length(prepped_data_EH.testspikes{i, i_cell}))+40, 'k.')
    end
    pause()
    clf
end