path = '/Volumes/Lab/Users/Nora/pillow_test/';
files = dir(path);

load('testmovie.mat');

datarun = load_data('2005-04-26-1/data004-from-data005-pca/data004-from-data005-pca');
datarun = load_neurons(datarun);

datarun_mas=load_data('2005-04-26-1/data005-pca/data005/data005');
datarun_mas=load_params(datarun_mas);

all_trial_starts = datarun.triggers(1:12:end);
trial_starts = all_trial_starts(1:25);

for i_file = 1:length(files)
    filename = files(i_file).name;
    if ~strcmp(filename(1),'.') && strcmp(filename(end), 't')
        cid = str2double(filename(1:(end-4)));
        load([path filename])
        xval = pillow_predict(fittedGLM, cid, testmovie_color, datarun, datarun_mas, trial_starts);
        plotraster(xval, fittedGLM, 'raster_length', 10);
        set(gcf,'PaperPositionMode','auto');
        set(gcf,'PaperOrientation','landscape');
        set(gcf, 'PaperSize', [20 4]);    
        set(gcf, 'Position', [0 0 1500 200])
        print(gcf,'-dpsc2',[path num2str(cid) '.eps']);
        close all
    end
end