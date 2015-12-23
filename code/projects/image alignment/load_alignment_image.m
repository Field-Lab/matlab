function alignment_image = load_alignment_image(ai_type,spec)
% load_alignment_image     load alignment image and x/ydata (its location in array or monitor coordinates)
%
% usage:  alignment_image = load_alignment_image(ai_type,spec)
%
% arguments:     ai_type - type of alignment image
%                           'pm'    - photographic mapping image
%
%                   spec - specifications for the image
%                           see the code!
%
% outputs:     alignment_image - struct with fields
%                       alignment_image.im    - matrix of the image
%                       alignment_image.xdata - horizontal coordinates of image (1st and last column)
%                       alignment_image.ydata - vertical coordinates of image (1st and last row)
%
%
%  NOTE: if called with no arguments, it returns a list of valid image names
%
%
% 2010-01  gauthier
%



% if called with no arguments, return valid image names
if nargin == 0
    alignment_image = {'pm32','pm10','pm2'};
    return
end


% identify which image to load
switch ai_type
    case 'pm'
        switch spec
            case 'pm32'  % BW white noise frame with 32 pixels/stixel
                image_path = [matlab_code_path 'code/projects/image alignment/images/photographic mapping 32.tif'];
                xdata = [1 640];  ydata = [81 400];
                
            case 'pm10'  % BW white noise frame with 10 pixels/stixel
                image_path = [matlab_code_path 'code/projects/image alignment/images/photographic mapping 10.tif'];
                xdata = [1 640];  ydata = [81 400];
                
            case 'pm2'  % BW white noise frame with 2 pixels/stixel
                image_path = [matlab_code_path 'code/projects/image alignment/images/photographic mapping 2.tif'];
                xdata = [1 640];  ydata = [81 400];
            otherwise
                error('alignment image specification ''%s'' not recognized.',num2str(spec))
        end
    otherwise
        error('alignment image type ''%s'' not recognized.',ai_type)
end

% save in output struct
alignment_image = struct;
alignment_image.im = imread(image_path);
alignment_image.xdata = xdata;
alignment_image.ydata = ydata;
