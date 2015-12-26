%This SCRIPT takes elecresp files and changes the paths so that analysis is possible outside of the
%usual filesystem structure of '/Volumes/Analysis/...'

%user chooses local folder with elecresp files in it. All .mat files in this folder that start with "elecresp" will be converted
dirn = uigetdir;
if dirn == 0; return; end
pdirn =	fileparts(dirn);
eidir = inputdlg('Enter EI folder data number', 'EI Location', 1, {'0'});
eidir = ['data' sprintf('%3.3d',str2num(eidir{1}))]; 
fileList = dir([dirn filesep '*mat']);

for i = 1:size(fileList, 1)

    clear elecResp;

    fn = [dirn filesep fileList(i).name];
    [~,fnname, fnext] = fileparts(fn);
    load(fn)

    %change the names
    elecResp.names.data_path = [dirn filesep];
    elecResp.names.rrs_ei_path = [pdirn filesep eidir filesep eidir '.ei'];
    elecResp.names.rrs_params_path = [pdirn filesep eidir filesep eidir '.params'];
    elecResp.names.savePath = [dirn filesep];
    save([dirn filesep fnname '_win' fnext], 'elecResp');
end
