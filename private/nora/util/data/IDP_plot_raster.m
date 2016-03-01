function  IDP_plot_raster(prepped_data, cell_to_plot)

% Visual check
figure;
hold on
for i = 1:size(prepped_data.testspikes, 1)
    plot(prepped_data.testspikes{i, cell_to_plot}, i*ones(length(prepped_data.testspikes{i, cell_to_plot})), 'k.')
end

end