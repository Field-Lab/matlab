function datarun = setsimplemarks(datarun, rgc, thresh)

sta = get_sta(datarun, rgc);
maxsta = max(mean(abs(sta), 3), [], 4);
simplemarks = maxsta > (thresh*robust_std(sta(:)));
datarun.stas.rfs{datarun.cell_nums(rgc)} = []; % Hack to make sure it's clear
datarun.stas.marks{datarun.cell_nums(rgc)} = simplemarks;