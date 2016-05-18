stim = {'WN', 'NSEM'};
files = dir(['/Volumes/Lab/Users/Nora/GLMFits/201208093/' stim{2} '/*Mid*.mat']);
Scores = zeros(2,length(files));
for i_cell = 1:length(files)
    for i_stim = [2 1]
        try
            load(['/Volumes/Lab/Users/Nora/GLMFits/201208093/' stim{i_stim} '/' files(i_cell).name])
            Scores(i_stim,i_cell) = fittedGLM.xval.corr;
        catch
            disp(['Cell ' files(i_cell).name ' does not have a ' stim{i_stim} ' fit'])
        end
    end
end
plot(Scores(1,:), Scores(2,:), '.k');
hold on; plot([0 1], [0 1]);
axis square