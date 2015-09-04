% # 2011-06-30-0 # big 4, high-res
% 
% Big 4
% 2005-05-26-7
% 2005-08-03-0
% 2007-03-27-1
% 2008-05-13-3
% 
% # Big 5 (OFF midgets are limiting)
% 2005-05-31-1
% 2013-05-28-4


path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/';
for_use = 'data014';

datarun = load_data(fullfile(path2data, for_use, for_use));
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta',[]);
datarun = get_sta_fits_from_vision(datarun);
datarun = load_sta(datarun,'load_sta', 'all');
datarun = set_polarities(datarun);

clear diam tc
for i = 1:4
    cellIDs = get_cell_indices(datarun, {i});
    
    for j=1:length(cellIDs)
        diam(j,i) = min(datarun.vision.sta_fits{j}.sd);
        if datarun.stas.polarities{j}<0
            tc(j,i) = min(datarun.vision.timecourses(j).g);
        else
            tc(j,i) = max(datarun.vision.timecourses(j).g);
        end
    end
end

figure
plot(tc(:),diam(:), '*')
axis([-0.3 0.3 0 10])

% cell 3002 (ON parasol) is a good candidate
