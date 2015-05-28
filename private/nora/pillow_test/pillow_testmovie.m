% % Load Movie
disp('Loading Stimulus Movies')
xml_file = '/Volumes/Analysis/stimuli/white-noise-xml/BW-20-1-0.48-11111.xml';
datarun = load_data('2005-04-26-1/data004-from-data005-pca/data004-from-data005-pca');
datarun = load_neurons(datarun);
testframes = 1140;

[temp_fitmovie,height,width,~,~] = get_movie(xml_file, datarun.triggers, testframes);
fitmovie_color=zeros(height,width,3,testframes);
for i=1:testframes
    fitmovie_color(:,:,:,i)=temp_fitmovie(:,:,:,ceil(i/interval));
end
disp(size(fitmovie_color))
clear temp_fitmovie height width itest