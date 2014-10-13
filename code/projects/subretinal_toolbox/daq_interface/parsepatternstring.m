function [M, param_names] = parsepatternstring(string)
% [M, param_names] = parsepatternstring(string)
% 
% This function parses a Pattern.ToString string into two cell arrays that
% can be used to process a dataset, or to load a stimulus from a logfile.
% It is spaghetti code, but dealing with string in Matlab is messy... 
% It could be worth switching to xml stimuli, instead of using strings to
% summarize strings.
%
% Parameters:
%   - string: the string to parse
%
% Returns:
%   - M: cell array with the relevant data parsed 
%   - param_names: name of each of the cells of M

% Version: 0.1 - 2012/09/18
% Author: Georges Goetz, Stanford University
% 

M = {};
param_names = {'Pattern type', 'Images data'};

delimiters = strfind(string, ';');

% Pattern type
M{1} = string(1:delimiters(1)-1);

% Parsing the image information
images_string = string(delimiters(1)+2:end);
delimiters = strfind(images_string, '}');
num_images = length(delimiters);
M{2} = [];
for kk=1:num_images
    if kk==1
        start_index = 2;
    else
        start_index = delimiters(kk-1) + 3;
    end
    end_index = delimiters(kk)-1;
    
    M{2} = [M{2} parseimagestring(images_string(start_index:end_index))];
end

    function image_struct = parseimagestring(string)
        % This function parses a sub-string representing a single image
        
        c_delimiters = strfind(string, ',');
        
        image_name_str = string(1:c_delimiters(1)-1);
        ind = find(image_name_str == '/'|image_name_str == '\', 1, 'last');
        image_struct.image = image_name_str(ind+1:end);
        
        image_struct.h_offset = str2double(string(c_delimiters(1)+11:c_delimiters(2)-1));
        if length(c_delimiters)==2
            image_struct.v_offset = str2double(string(c_delimiters(2)+11:end));
            image_struct.time = [];
        else
            image_struct.v_offset = str2double(string(c_delimiters(2)+11:c_delimiters(3)-1));
            image_struct.time = str2double(string(c_delimiters(3)+7:end-3));
        end
        
    end

end % parsepatternstring