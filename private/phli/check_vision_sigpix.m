datarun = load_data('2012-04-13-1/streamed/data001/data001');
datarun = load_params(datarun, 'keep_java_params', true);
datarun = load_neurons(datarun);
datarun = load_sta(datarun);

%%
datarun = get_sta_summaries(datarun, 'all', 'robust_std_method', 6);


%%
% Get Vision's version of pixel significance
for cid = datarun.cell_ids
    sta = datarun.stas.java_sta.getSTA(cid);
    mf = sta.getMainFrame();
    buf = permute(reshape(mf.getBuffer, 3, 60, 60), [3 2 1]);
    err = permute(reshape(mf.getErrorBuffer, 3, 60, 60), [3 2 1]);
    datarun.stas.sig{get_cell_indices(datarun, cid)} = sqrt(sum(buf.^2, 3) ./ sum(err.^2, 3));
end

% Confirm for CellID 124 (cell num 4) that the above calculation matches
% Vision's pixel significance calculation
sta = datarun.stas.java_sta.getSTA(124);
mf = sta.getMainFrame();
for j = 1:60
    for i = 1:60
        sig(j,i) = mf.getPixelSignificance(i-1,j-1);
    end
end
all(sig(:) == datarun.stas.sig{4}(:))

clear sta mf buf err sig i j cid ans


%% Plot comparisons
for cid = datarun.cell_types{4}.cell_ids
    figure;
    cnum = get_cell_indices(datarun, cid);
    
    sig = datarun.stas.sig{cnum};
    marks = datarun.stas.marks{cnum};
    
    sanesubplot(2, 3, {1,1});
    plot_rf(datarun, cid);
    
    sanesubplot(2, 3, {1,2});
    imagesc(marks);
    colormap gray;
    axis equal tight
    title 'matlab'
    
    sanesubplot(2, 3, {2,1});
    imagesc(sig > 3);
    colormap gray;
    axis equal tight
    title 3
    
    % Get the Vision significance threshold that will give the same number
    % of sigpix as the Matlab marks
    nmarks = sum(marks(:));
    nth_el = sig(:);
    nth_element_ip_parallel(nth_el, 60*60-nmarks);
    thresh = nth_el(60*60-nmarks);
    sigmatch = sig > thresh;
    
    sanesubplot(2, 3, {2,2});
    imagesc(sigmatch);
    colormap gray;
    axis equal tight
    title(num2str(thresh));
    
    sanesubplot(2, 3, {1,3});
    plot(datarun.stas.time_courses{cnum});
    
    sanesubplot(2, 3, {2,3});
    plot(datarun.vision.timecourses(cnum).g);
end


%%
% 49 305 558 651