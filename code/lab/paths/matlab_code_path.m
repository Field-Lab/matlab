function the_path = matlab_code_path
% matlab_code_path     return path to matlab code
%
% usage:  the_path = matlab_code_path()
%
%
% 2010-01  gauthier
% 2014-10  grosberg - transition to git file structure

the_path = fullfile(mfilename('fullpath'),'../../../../');