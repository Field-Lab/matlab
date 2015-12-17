function filter_matrix = showfitGMLM(fitGMLM2,text,mask,varargin)


p = inputParser;

% specify list of optional parameters
p.addParamValue('proj', eye(numel(mask(:))));
p.addParamValue('preproj', eye(numel(fitGMLM2.Linear.filter{1})));


% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

% assign variables from parsing.
proj = params.proj;
preproj = params.preproj;

nFilters = length(fitGMLM2.Linear.filter);

sta_dim1 = size(mask,1);
sta_dim2 = size(mask,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask));

filter_matrix = zeros(length(masked_frame),nFilters);

figure;
for ifilt=1:nFilters
subplot(2,2,ifilt)
xx = preproj*(fitGMLM2.Linear.filter{ifilt});
u_spatial = reshape_vector(xx(1:length(masked_frame)),masked_frame,indexedframe);
u_spatial = reshape(proj*u_spatial(:),[size(u_spatial,1),size(u_spatial,2)]);
imagesc((u_spatial));
%caxis([-0.3,0.3]);
colormap gray
colorbar
title(sprintf('%d',ifilt));

filter_matrix (:,ifilt) =xx(1:length(masked_frame));
end
suptitle(sprintf('%s',text));
end