% read data
data_path = '/Volumes/Palace/Data/gfield/2011-07-05-4/data008/';
data_file = 'data008000.bin';
samples = 2000;
channels = 520;
tic
data_matrix = read_vision_bin_file(data_path, data_file, samples, channels);
toc

figure(1)
plot(data_matrix(5,:))

% write data
file_name = 'test-data.bin';
tic
write_vision_bin_file(data_matrix, data_path, file_name);
toc

% read the written data to check it matches original
reloaded_data_matrix = read_vision_bin_file(data_path, file_name, samples, channels);

figure(2)
plot(reloaded_data_matrix(5,:))


%% test data with an odd number of channels

data_path = '/Volumes/Palace/Data/gfield/2009-04-13-5/data000/';
data_file = 'data000000.bin';
samples = 2000;
channels = 513;


data_matrix = read_vision_bin_file(data_path, data_file, samples, channels);
figure(1)
plot(data_matrix(3,:))


% write data
file_name = 'test-data.bin';
write_vision_bin_file(data_matrix, data_path, file_name);


% read the written data to check it matches original
reloaded_data_matrix = read_vision_bin_file(data_path, file_name, samples, channels);
figure(2)
plot(reloaded_data_matrix(11,:))


figure(3)
plot(reloaded_data_matrix(10,:))


%%


rdf = edu.ucsc.neurobiology.vision.io.RawDataFile('/Volumes/Palace/Data/gfield/2009-04-13-5/data000');
rdfh = rdf.getHeader();

%Get data:
   firstSample = 0;
   numSamples = 20000;
   data = rdf.getData(firstSample, numSamples);

   figure(2)
   plot(data(:,3))
   
   
%Write data back out:
   outpath = 'Users/peterli/MATLAB'
   % The data saver uses the header information and decides to save to
   % /Users/peterli/MATLAB/2011-07-14-6/data001.bin
   % I didn't mean for it to do this, but seems reasonable.

   edu.ucsc.neurobiology.vision.matlab.Matlab.saveRawData('Users/peterli/MATLAB', rdfh, data);

%Check that the data saved out look right:
   tdf = edu.ucsc.neurobiology.vision.io.RawDataFile('/Users/peterli/MATLAB/2011-07-14-6');
   all(tdf.getData(0,1) == rdf.getData(0,1))
   all(tdf.getData(9999,1) == rdf.getData(9999,1))
   % Both return 1, so the first and 10000th samples match
