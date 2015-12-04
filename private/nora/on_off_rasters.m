datarun_class = load_data('2015-10-06-0/data000-data015-norefit/data000-from-data000_data001_data002_data003_data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015/data000-from-data000_data001_data002_data003_data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015');
datarun_class = load_params(datarun_class);
datarun = load_data('2015-10-06-0/data000-data015-norefit/data002-from-data000_data001_data002_data003_data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015/data002-from-data000_data001_data002_data003_data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015');
datarun = load_neurons(datarun);

%%
datarun_class = load_data('2015-09-23-0/data004-data017-norefit/data006-from-data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015_data016_data017/data006-from-data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015_data016_data017');
datarun_class = load_params(datarun_class);
datarun = load_data('2015-09-23-0/data004-data017-norefit/data007-from-data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015_data016_data017/data007-from-data004_data005_data006_data007_data008_data009_data010_data011_data012_data013_data014_data015_data016_data017');
datarun = load_neurons(datarun);

%%
prepped_data = interleaved_data_prep(datarun, [3600*2 3600], 55, 'cell_spec','Off Parasol nc10', 'visual_check', 1, 'datarun_class', datarun_class);
count = 0;
for i = 1:30
    for j = 1:30
        a = a + length(prepped_data.testspikes{i,j});
        count = count +1;
    end
end
firing_rate = a/(count*30)

%%
figure;
for j = 1:2
    subplot(2,1,j)
    hold on
    for i = 1:size(prepped_data.testspikes, 1)
        plot(prepped_data.testspikes{i, j}, i*ones(length(prepped_data.testspikes{i, j})), 'k.')
    end
    xlim([0 3])
end