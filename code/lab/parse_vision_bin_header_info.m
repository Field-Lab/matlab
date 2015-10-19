function header_info = parse_vision_bin_header_info(fid)

% parse the header information
header_info = zeros(1,28);
for header_counter = 1:28
    
    if header_counter == 11
        header_info(header_counter)
        temp_info = fread(fid, 1, 'uint32', (header_info(10)-4), 'b');
        header_info(header_counter) = temp_info;
    
    else 
        temp_info = fread(fid, 1, 'uint32', 'b');        
        header_info(header_counter) = temp_info;
    
    end
end