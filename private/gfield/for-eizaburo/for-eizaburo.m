% plantain
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';

% blueberry
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-26-2/data001-s6369-s9551/data001-s6369-s9551';

% apple
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2010-03-05-2/data013/data013';


% load data
datarun = load_data(path_and_name{1,1});
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_neurons(datarun);

figure(1); clf;

for ct = 1:4

    cell_type = {ct};
    temp_cell_indices = get_cell_indices(datarun, cell_type);

    %------
    num_spike_list = zeros(length(temp_cell_indices),1);

    for cll = 1:length(temp_cell_indices)
        num_spike_list(cll) = length(datarun.spikes{temp_cell_indices(cll)});
    end

    mean_spike_rate = mean(num_spike_list) ./ (datarun.duration)
    
    % ---


    % plot the firing rate over time with 1 min bin

    bin_size = 60;
    num_bins = (datarun.duration - 6396) ./ bin_size;
    %num_bins = (datarun.duration) ./ bin_size;
    td_spike_rates = zeros(length(temp_cell_indices), num_bins);
    for cll = 1:length(temp_cell_indices)

        for bn = 1:num_bins
            temp_bin_begin = 6396 + bin_size * (bn-1);
            temp_bin_end = 6396 + (bin_size * bn);
     %       temp_bin_begin = bin_size * (bn-1);
     %       temp_bin_end = (bin_size * bn);

            temp_num = length(find(datarun.spikes{temp_cell_indices(cll)} >= temp_bin_begin & datarun.spikes{temp_cell_indices(cll)} < temp_bin_end));
            td_spike_rates(cll, bn) = temp_num ./ bin_size;
        end
    end

    std_spikes = std(td_spike_rates, 1);
    mean_spikes = mean(td_spike_rates,1);
    var_spikes = var(td_spike_rates,0,2);

    subplot(2,2,ct)
    %plot(std_spikes./mean_spikes)
    [hist_Vars, hist_bins] = hist(var_spikes);
    bar(hist_bins, hist_Vars)
    %plot(mean_spikes, 'k')
    hold on
    %plot(mean_spikes + std_spikes, 'r')
    %plot(mean_spikes - std_spikes, 'r')
    xlabel(datarun.cell_types{ct}.name)
    %ylabel('spike rate')
    %axis([1 52 0 1.5])
    hold off
    
    file_name = ['blueberry-',datarun.cell_types{ct}.name]
    spike_rates = td_spike_rates;
    save(file_name, 'spike_rates')
    
end