function [W,kernel_norms] = make_cone_weights_matrix(rf_size,kernel_spec,kernel_colors,varargin)
% make_cone_weights_matrix     Generate matrix in which each column is a cone kernel (reshaped to a vector)
%
% usage:  [W,kernel_norms] = make_cone_weights_matrix(rf_size,kernel_spec,kernel_colors,<params>)
%
% arguments:     rf_size - RF dimensions, first is y size, second is x size
%            kernel_spec - struct array giving the parameters of each cone as follows
%                  kernel_spec(kk).radius - radius of gaussian (in units of RF pixels)
%                                 .center - center of gaussian (in RF coordinates)
%                                 .type   - cone type, 'L', 'M', or 'S' (must match the field names of kernel_colors
%          kernel_colors - struct specifying RGB triplet for each cone type, e.g. kernel_colors.L = [.2 .8 .1]
%               <params> - struct or list of optional parameters (see below)
%
%
% outputs:             W - sparse matrix in which each column is an RGB cone kernel (reshaped to a vector)
%           kernel_norms - vector giving the L2 norm of each kernel
%
%
%
% optional params, their default values, and what they specify:
%
% verbose               true         	show output
% eff_radius_scale      2               size of the region in which to compute the kernel
%                                           in units of gaussian sigmas.  pixels outside this radius are set to 0.
%                                           the kernel is sparse, so this parameters determines the size (in bytes) of W
%
% 2009-09 gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', true);
p.addParamValue('eff_radius_scale',2);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;



% note number of kernels
num_kernels = length(kernel_spec);


if params.verbose
    T=text_waitbar('Generating W...');
    start_time = clock;
end

% initialize W
%W = sparse(prod(rf_size)*length(fieldnames(kernel_colors)),num_kernels);
W = sparse(prod(rf_size)*3,num_kernels); % the hardcoded 3 makes this an image

% fill in each cone

% for each location...
for kk=1:length(kernel_spec) 
%     
%     kk=4276
%     figure
%     tt=1
%     for cnt=1:3:12*3
%         kk=4998+cnt;
%     kk=5002
%     kk=5005
%     kk=5008
    % make BW kernel
    bw_kern = make_gaussian('dim',2,'x_size',rf_size(2),'y_size',rf_size(1),'normalize','sum',...
        'center_radius',kernel_spec(kk).radius,'center',kernel_spec(kk).center,...
        'effective_radius',ceil(params.eff_radius_scale*kernel_spec(kk).radius));
    
%     subplot(4,3,tt)
%     imagesc(bw_kern)
%     title([int2str(kk),',   ',num2str(kernel_spec(kk).center)]) 
%     tt=tt+1;
%     end 
    
    % reshape to a single column
    bw_kern = reshape(bw_kern,[],1);

    % get rgb values
    rgb = kernel_colors.(kernel_spec(kk).type);
    % multiply out
    rgb_kern = bw_kern * rgb;
    % reshape
    rgb_kern = reshape(rgb_kern,[],1);
    % normalize sum
    rgb_kern = rgb_kern/sum(rgb_kern);
    % put into W
    W(:,kk) = rgb_kern;

    % update waitbar
    if params.verbose && mod(kk,100)==0;
        T=text_waitbar(T,kk/length(kernel_spec));
    end
    
end



if params.verbose
    fprintf('\n   done (%0.1f seconds)\n',etime(clock,start_time));
end


% compute the normalization factor for each kernel
switch 2
    case 1 % all at once
        kernel_norms = diag(W'*W)';
    case 2 % in chunks
        kernel_norms = zeros(1,num_kernels);
        chunk_size = 1e4;
        T=text_waitbar('Computing norm factors...');
        num_chunks = ceil(num_kernels/chunk_size);
        for cc = 1:num_chunks
            subset_start = 1 + (cc-1)*chunk_size;
            subset_end = min(subset_start + chunk_size - 1 , num_kernels);
            kernel_norms(subset_start:subset_end) = diag(W(:,subset_start:subset_end)'*W(:,subset_start:subset_end))';
            T=text_waitbar(T,cc/num_chunks);
        end
end 






% % make BW kernel
% bw_kern = make_gaussian('dim',2,'x_size',rf_size(2),'y_size',rf_size(1),'center_radius',kernel_radius,...
%     'center',kernel_centers(kk,:),'effective_radius',params.eff_radius);
% % reshape to a single column
% bw_kern = reshape(bw_kern,[],1);
%
% % generate this kernel in each color
% for cc = 1:kern_c
%     % get rgb values
%     rgb = kernel_colors.(kernel_color_names{cc});
%     % multiply out
%     rgb_kern = bw_kern * rgb;
%     % reshape
%     rgb_kern = reshape(rgb_kern,[],1);
%     % put into W
%     W(:,(kk-1)*3 + cc) = rgb_kern;
%
%     % plot each kernel
%     %figure(1);clf;imagesc(norm_image(reshape(rgb_kern,31,41,3)>0));pause
% end