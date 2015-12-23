function imagepath = get_array_image(array_id)
% GET_ARRAY_IMAGE
% usage: image_file = get_array_image(array_id)
%
% 2010-10 phli, abstracted out of load_array_info
%

if array_id < 500  % 64 electrode array
    image_file = '61 map.tif';
elseif (array_id >= 500) && (array_id < 1500) % 512 array
    image_file = '512 map.tiff';
elseif (array_id >= 1500) && (array_id < 2500) % 519 array
    image_file = '519 map.tiff';
else
    error('array id %d not recognized.', array_id)
end
imagepath = fullfile( 'code/projects/image alignment/images' ,image_file);