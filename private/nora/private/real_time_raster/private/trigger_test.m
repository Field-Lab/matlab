rawdatafiledir = '/Volumes/Stream/Data/9999-99-99-9/data008.bin';
rdf=edu.ucsc.neurobiology.vision.io.RawDataFile(rawdatafiledir);
trigger_data=rdf.getData(0,0,100*20000);
trigger_times = find(diff(trigger_data)>500);
hist(diff(trigger_times)/20000, 100)