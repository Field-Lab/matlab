threshold = 0.5;
n_blocks = 37;
for i_cell = 1:15
    res_spikes = [];
    count_all = zeros(1,50);
    for i_block = 1:n_blocks
        res_spikes_temp = lin_d{i_cell, i_block}(logical(xval{i_block}.rasters.recorded) & (xval{i_block}.glm_rateperbin < threshold*mean(xval{i_block}.glm_rateperbin)));
        res_spikes = [res_spikes res_spikes_temp];
        count_all_temp = hist(lin_d{i_cell,i_block}, center_all);
        count_all = count_all + count_all_temp;
    end
    count_res = hist(res_spikes, center_all);
    plot(center_all,count_all/(72000*n_blocks))
    hold on
    plot(center_all, count_res/length(res_spikes))
    legend('Drive at All Times', 'Drive at Residual Spike Times');
    pause()
    hold off
end

