function rf = rf_from_sta(sta, varargin)
% rf_from_sta compute spatial summary (RF) of an STA movie
%
% This function identifies the STA's overall timecourse, and computes the amplitude of
% this timecourse at each stixel. Both ON and OFF cells have positive centers
% and negative surrounds. For cells with strong ON and OFF components,
% (e.g. SBCs), the polarity will depend on user parameters.
% For a color STA, color can be in the timecourse (in which case the RF will be one dimensional),
% or color can be in the RF (in which case the timecourse will be one dimensional).
%
%
% usage: rf = rf_from_sta(sta, varargin)
%
% arguments: sta - standard 4D matrix of one STA (y,x,color,time)
% varargin - struct or list of optional parameters (see below)
%
% outputs: RF - STA spatial summary frame, matrix dimensions of height, width, and color
%
%
% optional fields in params, their default values, and what they specify:
%
% method 'project' how to get the RF
% 'project' - project the timecourse at each stixel
% onto the overall timecourse
% 'mean' - take the mean of a few frames
% 'std' - take the std of a few frames
% 'svd' - perform svd - NOT SUPPORTED YET
%
%
% if method == 'project'
% frames ':' which frames to use (see parse_frame_spec)
% sig_stixels [] matrix of significant stixels.
% If empty, significant_stixels is called with default parameters.
% num_colors size(sta,3) how many colors to have in the RF
% template 'max' what kind of overall timecourse to use
% interpretation depends on the value of num_colors:
% if num_colors == 1
% x project color x onto the overall timecourse from color x
% [x y z] Concatenate the overall timecourses from colors x, y, & z, and use them to
% project the concatenated timecourses of colors 1, 2, and 3
% this vector must be the same length as num_colors
%
% if num_colors == size(sta,3)
% generate a single timecourse, and project each color channel separately
% 'max' - get the overall timecourse from each color, and use the one with the largest variance
% 'sum' - get the overall timecourse from each color, and use the sum
% [x y z] - project color 1 onto overall timecourse x
% project color 2 onto overall timecourse y, etc.
% this vector must be the same length as num_colors
%
%
% if method == 'mean'
% frames ':' which frames to use (see parse_frame_spec)
% colors 1:size(sta,3) which colors to use
%
%
% if method == 'std'
% frames ':' which frames to use (see parse_frame_spec)
% colors 1:size(sta,3) which colors to use
%
%
% if method == 'svd'
% frames ':' which frames to use (see parse_frame_spec)
% colors 1:size(sta,3) which colors to use
% sig_stixels [] matrix of significant stixels.
% If empty, significant_stixels is called with default parameters.
%
%
%
%
% 2008-10 gauthier
% 2010-03 GDF: modified to allow function to handle STAs that are 1D in space
%
% SET UP OPTIONAL ARGUMENTS
p = inputParser;
% specify list of optional parameters
p.addParamValue('method', 'project', @(x)any(strcmpi(x,{'project','mean','std','svd'})));
% forked parameters
% method == 'project'
p.addParamValue('frames','default value');
p.addParamValue('sig_stixels','default value');
p.addParamValue('num_colors','default value');
p.addParamValue('template', 'default value');
% method == 'mean','std'
p.addParamValue('colors','default value');
% resolve user input and default values
p.parse(varargin{:});
% get params struct
params = p.Results;
% set parameters forked from method
% set default values in an inputParser object
p_temp = inputParser;
switch params.method
case 'project'
p_temp.addParamValue('frames',':');
p_temp.addParamValue('sig_stixels',[]);
p_temp.addParamValue('num_colors',size(sta,3));
% base the default 'template' value on the number of colors in the RF
if size(sta,3) > 1
p_temp.addParamValue('template', 'max');
else
p_temp.addParamValue('template', 1);
end
case 'mean'
p_temp.addParamValue('frames',':');
p_temp.addParamValue('colors',1:size(sta,3));
case 'std'
p_temp.addParamValue('frames',':');
p_temp.addParamValue('colors',1:size(sta,3));
case 'svd'
p_temp.addParamValue('sig_stixels',[]);
p_temp.addParamValue('colors',1:size(sta,3));
p_temp.addParamValue('frames',':');
end
% add forked parameters to the params struct
params = add_forked_params(params,p_temp);
% BODY OF THE FUNCTION
% if STA is empty, return an empty RF
if isempty(sta) || (~isempty(params.sig_stixels) && all(reshape(params.sig_stixels,[],1) == 0))
rf = [];
return
end
% convert to double to be sure the result is double precision
sta = double(sta);
% compute the RF
switch params.method
case 'project'
% parse frame specification
frames = parse_frame_spec(params.frames,size(sta,4));
% keep only desired frames
sta = sta(:,:,:,frames);
% get significant stixels, if needed
if isempty(params.sig_stixels)
params.sig_stixels = significant_stixels(sta);
end
% get timecourse of each color
timecourse = time_course_from_sta(sta,params.sig_stixels);
% if timecourse is empty, return an empty RF
if isempty(timecourse)
rf = [];
return
end
% zero mean the time course for each color independently
for num_clr = 1:length(timecourse(1,:))
timecourse(:,num_clr) = timecourse(:,num_clr) - repmat(mean(timecourse(:,num_clr)), length(timecourse(:,num_clr)), 1);
end
% zero mean the STA color by color !!! GDF modified this bit of code 2010-03-15
[nrw, ncl, nclr, nfrms] = size(sta);
for clr = 1:nclr
% calculate mean across frames and subtract this mean from color
if ndims(mean(squeeze(sta(:,:,clr,:)))) == 3
mean_clr = mean(squeeze(sta(:,:,clr,:)), 3);
sta(:,:,clr,:) = squeeze(sta(:,:,clr,:)) - repmat(mean_clr, [1 1,nfrms]);
elseif ndims(mean(squeeze(sta(:,:,clr,:)))) == 2 % handles case when STA is 1 dimensional in space
mean_clr = mean(squeeze(sta(:,:,clr,:)), 2);
sta(:,:,clr,:) = squeeze(sta(:,:,clr,:)) - repmat(mean_clr, [1, nfrms]);
end
end
% do the projection
switch params.num_colors
case 1 % one color in the RF
switch length(params.template)
case 1
% get STA values from a single color
sta_for_projection = reshape(sta(:,:,params.template,:),[],size(sta,4));
% get timecourse from the same color
tc_for_projection = timecourse(:,params.template)';
case size(sta,3)
% concatenate STA values from all colors
sta_for_projection = reshape(permute(sta,[1 2 4 3]),[],size(sta,3)*size(sta,4));
% concatenate timecourses from the desired colors
tc_for_projection = reshape(timecourse(:,params.template),[],1)';
otherwise
error('params.template must have length 1 or length %d.',size(sta,3))
end
rf = get_projection(sta_for_projection,tc_for_projection);
% reshape
rf = reshape(rf,size(sta,1),size(sta,2));
otherwise % multiple colors in the RF
switch class(params.template)
case 'char'
switch params.template
case 'sum' % project each color onto the sum of overall timecourses
% sum the timecourses
tc_for_projection = sum(timecourse,2)';
% reshape sta
sta_for_projection = reshape(sta,[],size(sta,4));
% compute projection
rf = get_projection(sta_for_projection,tc_for_projection);
% reshape
rf = reshape(rf,size(sta,1),size(sta,2),size(sta,3));
case 'max' % project each color onto the biggest variance timecourse
% get the variance of each timecourse
tc_vars = var(timecourse);
% identify largest
[junk,tc_max] = max(tc_vars);
% use it
tc_for_projection = timecourse(:,tc_max)';
% reshape sta
sta_for_projection = reshape(sta,[],size(sta,4));
% compute projection
rf = get_projection(sta_for_projection,tc_for_projection);
% reshape
rf = reshape(rf,size(sta,1),size(sta,2),size(sta,3));
otherwise
error('template type ''%s'' not recognized',params.template)
end
case 'double' % project specified colors onto the specified overall timecourses
% error checking
% ensure enough colors are specified
if length(params.template) ~= params.num_colors
error('params.template must have length params.num_colors == %d',params.num_colors)
end
% normalize vector length of each timecouse
timecourse = normalize_to_unit_sphere(timecourse')';
% initialize rf
rf = zeros(size(sta,1)*size(sta,2),params.num_colors);
% for each of the desired output colors
for cc=1:params.num_colors
% get the relevant timecourse
tc_for_projection = timecourse(:,params.template(cc))';
% get the relevant color
sta_for_projection = reshape(sta(:,:,cc,:),[],size(sta,4));
% compute projection
rf(:,cc) = get_projection(sta_for_projection,tc_for_projection);
end
% reshape
rf = reshape(rf,size(sta,1),size(sta,2),params.num_colors);
otherwise
fprintf('\n\n\n')
disp(params.template)
error('above template type not recognized')
end
end
case 'svd'
error('svd not implemented yet!')
case 'mean'
% get frames and color
% parse frame specification
frames = parse_frame_spec(params.frames,size(sta,4));
% average frames
rf = mean(sta(:,:,params.colors,frames),4);
case 'std'
% get frames and color
% parse frame specification
frames = parse_frame_spec(params.frames,size(sta,4));
% average frames
rf = std(sta(:,:,params.colors,frames),[],4);
otherwise
error('Method for getting RF ''%s'' not recognized.',params.method)
end
function dp = get_projection(matrix,vect)
dp = matrix*vect';