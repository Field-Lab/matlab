%%
axon = 9;
somdenfn = ['somden' num2str(axon)];
somdenpath = fullfile(server_path(), '/2007-09-18-4/somden', [somdenfn '.tif']);
somdenpath_despeck = fullfile(server_path(), '/2007-09-18-4/somden', [somdenfn '_despeck.tif']);

if exist(somdenpath_despeck, 'file')
    somdenpath = somdenpath_despeck;
end

%%
inf = imfinfo(somdenpath);
slices = regexp(inf(1).ImageDescription, 'slices=(\d+)', 'tokens', 'once');
slices = str2double(slices{1});

%%
stack = zeros(inf(1).Height, inf(1).Width, 3, slices);
for i = 1:slices
    % Assume G and B channels and composite
    stack(:,:,2,i) = imread(somdenpath, i*2 - 1);
    stack(:,:,3,i) = imread(somdenpath, i*2);
end

stack_gfp = squeeze(stack(:,:,2,:));
stack_gfp_sm = smooth3(stack_gfp, 'gaussian', [5 5 5]);


%%
zlevels{9} = [1:4:slices];
contours{9} = 40*[1 1];
zlevels{17} = [1:2:slices];
contours{17} = 25*[1 1];

cla; 
p = contourslice(stack_gfp_sm, [], [], zlevels{axon}, contours{axon}); 
set(gca, 'DataAspectRatio', 1./[.232 .232 .279], 'Color', 'none'); 
set(gcf, 'Color', 'k');
set(p, 'CDataMapping', 'direct');


%%
small = 15;
smallp = arrayfun(@(P) size(get(P, 'Vertices'), 1), p) < small; 
delete(p(smallp)); 
p = p(~smallp);

%%
zdata = arrayfun(@(P) first(get(P, 'ZData')), p);
ncolors = size(colormap,1);

% Scale to fill color space
colors = zdata;
colors = colors - min(colors);
colors = round(colors .* ncolors ./ max(colors));

for i = 1:length(zdata)
    set(p(i), 'CData', colors(i)*ones(length(get(p(i), 'Faces')), 1));
end