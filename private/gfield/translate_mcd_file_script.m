clear


%% get a header from a generic file: This will be used as the header for
%% the final file -- it will be an incorrect header
visvolatile();
filepath = '/Data/MCS/gfield/2009-04-13-5/data000/data000000.bin/';
rdf = edu.ucsc.neurobiology.vision.io.RawDataFile(filepath);
rdfh = rdf.getHeader();
rdfh.methods();

%%

% data reading path
cd /Data/MCS/gfield/2011-08-04-1/data007/

% data writing path
path_to_write = 'Volumes/Boxer/Data/gfield/2011-08-04-1/data007-80Hz/data007/';

% provide information for reading data
num_files = 25;
data_run_number = '007';
filename = cell(num_files,1);
for cnt = 1:num_files;
    if cnt-1 < 10
        filename{cnt} = ['data',data_run_number, '00',num2str(cnt-1),'.mcd'];
    elseif cnt-1 >=10 && cnt-1 < 100
        filename{cnt} = ['data',data_run_number, '0',num2str(cnt-1),'.mcd'];
    else
        filename{cnt} = ['data',data_run_number,num2str(cnt-1),'.mcd'];
    end
end



path_to_el_correspondences = '~/MCS-conversion/mapping.mat';
load(path_to_el_correspondences)

    

%%

% set library path
[pathname, name, ext]=fileparts(which('nsMCDLibrary.dylib'));
ns_SetLibrary([pathname filesep name ext]);


% Open data file
[~, hfile] = ns_OpenFile(filename{1});

[~, FileInfo] = ns_GetFileInfo(hfile);


% get sampling rate from the info on the first electrode
[nsresult, analog_file_info] = ns_GetAnalogInfo(hfile, 1);
sampling_rate = analog_file_info.SampleRate;

% setup filters
% highpass
filter_order = 6;
filter_cutoff = 80; % Hz
fNorm = filter_cutoff ./ (sampling_rate/2);
[b,a] = butter(filter_order, fNorm, 'high');

% lowpass
filter_order = 6;
filter_cutoff = 5000; % Hz
fNorm = filter_cutoff ./ (sampling_rate/2);
[b_low,a_low] = butter(filter_order, fNorm, 'low');


num_channels = FileInfo.EntityCount;
if num_channels == 505
    channel_list = [1:252,505];
elseif num_channels == 253
    channel_list = 1:1:253;  % (503) unfiltered data
else
    error('entry in FileInfo.EntityCount is not recognized, check data file')
end
ns_result = ns_CloseFile(hfile);

% Fudge some parameters
%file = 'Data/gfield/2011-08-04-1/data000/'; % no starting slash

commonpath = ''; % Not needed as we are only saving to one location
samples_to_buffer = 100; % Somewhat arbitrary but this seems to work okay
num_buffers = 250; % Somewhat arbitrary but this seems to work okay
seconds_to_stream = 0; % Don't stop
% Open RawDataSaver
rds = edu.ucsc.neurobiology.vision.io.RawDataSaver(path_to_write, commonpath, rdfh, samples_to_buffer, num_buffers, seconds_to_stream);


for chunk = 1:length(filename);
    
    fprintf('processing chunk %d of %d \n',chunk, length(filename)) 
 
    % Open data file
    [~, hfile] = ns_OpenFile(filename{chunk});

    [nsresult, FileInfo] = ns_GetFileInfo(hfile);
    record_duration = FileInfo.TimeSpan * sampling_rate;
    
    piece_size = 400000;
    num_pieces = floor(record_duration ./ piece_size) + 1;
    remaining = mod(record_duration, piece_size);

    first_half = floor(record_duration ./ 2);
    second_half = first_half+1;

    
    for cnt = 1:num_pieces
        if cnt == num_pieces
            output_matrix = zeros(remaining, 513);
            [nsresult, ContinuousCount, output_matrix(:,1:253)] = ns_GetAnalogData(hfile, channel_list, (piece_size*(cnt-1))+1, remaining);
            
        else
            output_matrix = zeros(piece_size, 513);
            [nsresult, ContinuousCount, output_matrix(:,1:253)] = ns_GetAnalogData(hfile, channel_list, (piece_size*(cnt-1))+1, piece_size);
        end
            
        % rescale the data
        output_matrix(:,1:252) = output_matrix(:,1:252) .* 10000000;

        % move TTL signal to first row of matrix
        wave_indices = [253, 1:252];
        output_matrix(:,1:253) = output_matrix(:,wave_indices);

        %rescale TTL signal
        TTL_sig = output_matrix(:,1);
        low_indices = find(TTL_sig < 0);
        high_indices = find(TTL_sig > 3);
        TTL_sig(low_indices) = 0;
        TTL_sig(high_indices) = -2048;
        ttl_inds = find(TTL_sig == -2048);
        ii = 1;
        while ii < length(ttl_inds)
            if TTL_sig(ttl_inds(ii)) == (TTL_sig(ttl_inds(ii)-1) -2048)

                TTL_sig(ttl_inds(ii)+1) = -2048;
                TTL_sig(ttl_inds(ii)+2) = -2048;
                TTL_sig(ttl_inds(ii)+3) = 0;
                ii = ii + 4;
            else
                ii = ii +1;
            end
        end
        output_matrix(:,1) = TTL_sig;

        % put in noise on filler channels
        %output_matrix(:,254:513) = normrnd(0, 10, record_duration, 260);
        if cnt == num_pieces
            output_matrix(:,254:513) = normrnd(0, 10, remaining, 260);
        else
            output_matrix(:,254:513) = normrnd(0, 10, piece_size, 260);
        end
      
            
        % filter the data with the highpass filter
        output_matrix(:,2:513) = filtfilt(b_low, a_low, output_matrix(:,2:513));
        output_matrix(:,2:513) = filtfilt(b, a, output_matrix(:,2:513));


        % force channels to integers
        output_matrix = round(output_matrix);

        % rearrange channels to correspond to rrs electrode map
        new_output_matrix = zeros(size(output_matrix));
        new_output_matrix(:,mapping) = output_matrix;


        % Use RawDataSaver to save iteratively
        edu.ucsc.neurobiology.vision.matlab.Matlab.saveRawData(rds, new_output_matrix); % samples 0-999
    end
    ns_Result = ns_CloseFile(hfile);

end

% Close RawDataSaver; important
rds.finishSampleProcessing();


