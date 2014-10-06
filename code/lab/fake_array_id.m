function array_id = fake_array_id(array_type)
% FAKE_ARRAY_ID    Return a valid array idea for a given array type
% usage: arr_id = fake_array_id(array_type)
%
% There are a good number of situations where it doesn't matter which array
% id you use as long as it's the right array type (i.e. 61, 512, or 519 
% electrode).  This will make sure you get a valid array_id for a given
% array_type.
%
% 2010-05 phli
%

switch array_type
    case 61
        array_id = 1;
    case 512
        array_id = 504;
    case 519
        array_id = 1501;
    otherwise
        array_id = [];
end