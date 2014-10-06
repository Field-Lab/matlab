function data_matrix = read_vision_bin_file(data_path, data_file, samples, channels)
%
% function data_matrix = read_vision_bin_file('data_path', samples)
%
% Reads a .bin file into MATLAB. 
%
% Inputs
%   data_path   character string of file path to data
%   data_file   char string name of file to read into MATLAB
%   samples     the number of samples to read in
%   channels    the number of channels in the recording
%
% Outputs
%   data_matrix     MxN matrix with M = channel number and N = sample
%                   number
%
% Author: GDF
% Date: 2011-07-20
%

cd(data_path)

%file_info = dir(data_file);
%total_bytes = file_info.bytes;

fid = fopen(data_file, 'r', 'b');
header_info = parse_vision_bin_header_info(fid);


% header is sometimes larger depending on when the data file was generated,
% so compare header size as indicated in header to this size and scan
% forward until the file pointer is at data
header_size = header_info(3);
if header_size > ftell(fid)
    byte_difference = header_size - ftell(fid);
    fseek(fid, byte_difference, 0);
end

% header stores number of channels, check to make sure this equals that
% specified by user
header_channel_number = header_info(17);
if channels ~= header_channel_number
    disp('Warning: number of channels ~= that in file header, using user value')
end
channel_number = channels;




data_matrix = zeros(channel_number,samples);

disp('reading data...');


if mod(channel_number,2) == 0 % if channel # is even    

    odd_channels = 1:2:channel_number;
    even_channels = 2:2:channel_number;
    byte_number = channel_number * samples * 3/2;

    byte_stream = fread(fid, byte_number, 'uint8','b');
  
    byte_stream = reshape(byte_stream, 3, []);

    b1 = reshape(byte_stream(1,:), channels/2, samples);
    b2 = reshape(byte_stream(2,:), channels/2, samples);
    b3 = reshape(byte_stream(3,:), channels/2, samples);

    % decompress unsigned integers into channels
    odd_channel_data = bitshift(b1, 4) + bitshift(b2, -4) - 2048;
    even_channel_data = bitshift(bitand(b2, hex2dec('000F')), 8)  + b3 - 2048;

    data_matrix(odd_channels,:) = odd_channel_data;
    data_matrix(even_channels,:) = even_channel_data;

    
else %  if channel # is odd
   
    TTL_channel = samples * 2;
    odd_channels = 3:2:channel_number;
    even_channels = 2:2:channel_number;
    byte_number = TTL_channel + ((channel_number-1) * samples * 3/2);

    byte_stream = fread(fid, byte_number, 'uint8','b');
  
    % extract TTL signals
    % indices to TTLs
    TTL_indices_one = 1:ceil(channel_number*3/2):byte_number;
    TTL_indices_two = 2:ceil(channel_number*3/2):byte_number;
    TTL_indices =[TTL_indices_one; TTL_indices_two];
    TTL_indices = reshape(TTL_indices, 1,[]);

    TTL_bytes = byte_stream(TTL_indices);
    
%     TTL_bytes(1:100)
%     pause
%     X = find(TTL_bytes ~= 0)
%     TTL_bytes(X)
%     pause
    
    TTL_signals = typecast(uint8(TTL_bytes), 'uint16');
    
    
    % remove TTL signals from byte_stream
    byte_indices = true(byte_number,1);
    byte_indices(TTL_indices) = false;
    byte_stream = byte_stream(byte_indices);

    % reshape to extract each byte-triplet
    byte_stream = reshape(byte_stream, 3, []);

    b1 = reshape(byte_stream(1,:), (channels-1)/2, samples);
    b2 = reshape(byte_stream(2,:), (channels-1)/2, samples);
    b3 = reshape(byte_stream(3,:), (channels-1)/2, samples);

    % decompress unsigned integers into channels
    even_channel_data = bitshift(b1, 4) + bitshift(b2, -4) - 2048;
    odd_channel_data = bitshift(bitand(b2, hex2dec('000F')), 8)  + b3 - 2048;
    
    % fill in the data matrix with TTL and electrode data
    data_matrix(1,:) = TTL_signals;
    data_matrix(odd_channels,:) = odd_channel_data;
    data_matrix(even_channels,:) = even_channel_data;
    
end

fclose(fid);