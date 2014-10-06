function write_vision_bin_file(data_matrix, data_path, file_name, varargin)
%
% write_vision_bin_file(data_matrix, data_path, file_name, varargin)
%
% Give the function a data matrix, data_path, and file_name and it will
% write the data to a .bin file readable by vision.  The data must be
% 16-bit signed integers.  It 
%
%
% inputs
%   data_matrix   MxN matrix      M channels, N samples
%   data_path     char            character string of data path for saving
%                                   bin file
%   file_name     char            name of file to write
%
% Optional inputs for header
%   write_header        true        logical for whether to write a header for
%                                       file
%   array_id            1501        id for electrode array
%   trigger_delay       1000        delay before trigger (s * 1000000)
%   trigger_interval    75000       trigger interval (s * 1000000)
%   duration            900         duration of recording (s)
%   append_flag         false       false = make new file, true = append to
%                                       file
%
% Author GDF
% date: 2011-07-20

p = inputParser;

p.addParamValue('sampling_interval', 20000, @isnumeric)
p.addParamValue('array_id', 1501, @isnumeric);
p.addParamValue('trigger_delay', 1000, @isnumeric);  % s * 1000000
p.addParamValue('trigger_interval', 75000, @isnumeric); % s * 1000000
p.addParamValue('duration', 900, @isnumeric); % duration of experiment (s)
p.addParamValue('append_flag', false, @islogical); % whether the append to end of file

p.parse(varargin{:})

sampling_interval = p.Results.sampling_interval;
array_id = p.Results.array_id;
trigger_delay = p.Results.trigger_delay;
trigger_interval = p.Results.trigger_interval;
duration = p.Results.duration;
append_flag = p.Results.append_flag;

cd(data_path);

% get the number of channels and samples
[num_channels, num_samples] = size(data_matrix);


%% initialize and open file
% open file for writing (w), big-endian format (b)
if append_flag
    fid = fopen(file_name, 'a', 'b');
else
    fid = fopen(file_name, 'w', 'b'); 
end


%% write header?
if ~append_flag
    % header tag
    fwrite(fid, 0, 'uint32', 'b'); % tag = 0
    fwrite(fid, 4, 'uint32', 'b'); % 4 bytes
    fwrite(fid, 184, 'uint32', 'b'); % header length in bytes
    
    % time tag
    fwrite(fid, 1, 'uint32', 'b'); % tag = 1
    fwrite(fid, 12, 'uint32', 'b'); % 12 bytes
    fwrite(fid, 1904, 'uint32', 'b'); % 1904 = base year
    fwrite(fid, 0, 'uint32', 'b');
    fwrite(fid, 3e9, 'uint32', 'b');
    
    % user header
    fwrite(fid, 2, 'uint32', 'b'); % tag = 2
    fwrite(fid, 48, 'uint32', 'b'); % size of header in bytes, then skip 48 bytes
    fwrite(fid, 0, 'uint32', 44, 'b'); % size of header in bytes, then skip 48 bytes
    
    % format
    fwrite(fid, 3, 'uint32', 'b'); % tag = 3
    fwrite(fid, 4, 'uint32', 'b'); % 4 bytes
    fwrite(fid, 1, 'uint32', 'b'); % format = 1 for 25% compression
    
    % electrode array
    fwrite(fid, 4, 'uint32', 'b'); % tag = 4
    fwrite(fid, 8, 'uint32', 'b'); % 8 bytes
    fwrite(fid, num_channels, 'uint32','b'); % channel number
    fwrite(fid, array_id, 'uint32', 'b');
    
    % sampling frequency
    fwrite(fid, 5, 'uint32', 'b'); % tag = 5
    fwrite(fid, 4, 'uint32', 'b'); % 4 bytes
    fwrite(fid, sampling_interval, 'uint32', 'b'); % sampling frequency
    
    % trigger parameters
    fwrite(fid, 6, 'uint32', 'b'); % tag = 6
    fwrite(fid, 8, 'uint32', 'b'); % 8 bytes
    fwrite(fid, trigger_delay, 'uint32', 'b'); % first trigger delay in s * 1000000
    fwrite(fid, trigger_interval, 'uint32', 'b'); % first trigger interval in s * 1000000
    
    % header end
    fwrite(fid, 499, 'uint32', 'b'); % tag = 499
    fwrite(fid, 4, 'uint32', 'b'); % 4 bytes
    fwrite(fid, duration, 'uint32', 'b'); % number of seconds of recording
    
    % add some buffer to the header for extra space
    fwrite(fid, 0, 'uint32', 24,'b');
end

%% Write out data w/ compression



if mod(num_channels,2) == 0 % for an even number of channels

    data_matrix = uint16(data_matrix + 2048);

    % even and odd channel indices
    odd_channels = 1:2:num_channels;
    even_channels = 2:2:num_channels;

    % convert 2 12-bit numbers into 3 8-bit numbers
    b1 = bitshift(data_matrix(odd_channels,:), -4);
    b2 = bitshift(bitand(data_matrix(odd_channels, :), hex2dec('000F')), 4) + bitshift(data_matrix(even_channels, :), -8);
    b3 = bitand(data_matrix(even_channels, :), hex2dec('00FF'));

    % reshape matrices into vectors
    b1 = reshape(b1, 1, []);
    b2 = reshape(b2, 1, []);
    b3 = reshape(b3, 1, []);

    % combine the 3 8-bit components
    byte_stream = [b1; b2; b3];

    % reshape again into a vector
    byte_stream = reshape(byte_stream, 1, []);

    % write data
    fwrite(fid, byte_stream, 'uint8', 'b');

else % for an odd number of channels   

    data_matrix_no_TTL = uint16(data_matrix(2:num_channels,:) + 2048);
  
  
    % initialize the byte_stream
    byte_number = (2*num_samples) + ((num_channels-1) * num_samples * 3/2);
    byte_stream = zeros(byte_number, 1,'uint8');
 
    % channel indices
    odd_channels = 1:2:(num_channels-1);
    even_channels = 2:2:(num_channels-1);
    
    % TTL _data
    TTL_data = data_matrix(1,:);
    TTL_byte_data = typecast(uint16(TTL_data), 'uint8');
if 0    
    find(TTL_data ~= 0)
    
    length(TTL_byte_data)
    find(TTL_byte_data ~= 0)
    
    pause
end   
    %TTL indices
    TTL_indices_one = 1:((num_channels-1)*3/2)+2:byte_number;
    TTL_indices_two = 2:((num_channels-1)*3/2)+2:byte_number;
    TTL_indices =[TTL_indices_one; TTL_indices_two];
    TTL_indices = reshape(TTL_indices, 1,[]);    

    
    % convert 2 12-bit numbers into 3 8-bit numbers
    b1 = bitshift(data_matrix_no_TTL(odd_channels,:), -4);
    b2 = bitshift(bitand(data_matrix_no_TTL(odd_channels, :), hex2dec('000F')), 4) + bitshift(data_matrix_no_TTL(even_channels, :), -8);
    b3 = bitand(data_matrix_no_TTL(even_channels, :), hex2dec('00FF'));
    
    % reshape matrices into vectors
    b1 = reshape(b1, 1, []);
    b2 = reshape(b2, 1, []);
    b3 = reshape(b3, 1, []);

    % combine the 3 8-bit components
    electrode_stream = [b1; b2; b3];

    % reshape again into a vector
    electrode_stream = reshape(electrode_stream, 1, []);

    % get logical indices for non-TTL channels 
    byte_indices = true(1,byte_number);
    byte_indices(TTL_indices) = false;    
    
    
    % put TTL data and electrode data into byte_stream
    byte_stream(TTL_indices) = TTL_byte_data;
    byte_stream(byte_indices) = electrode_stream;
    
    
    % write data
    fwrite(fid, byte_stream, 'uint8', 'b');
    
end
    
% close file
fclose(fid);

   



