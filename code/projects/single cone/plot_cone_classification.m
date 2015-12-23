function plot_cone_classification(datarun, varargin)
% plot_cone_classification     Make cone classification plots
%
% usage:  plot_cone_classification(datarun, varargin)
%
% arguments:     datarun - datarun struct
%               varargin - struct or list of optional parameters (see below)
%
% outputs:     result - result of computation
%
%
% optional params, their default values, and what they specify:
%
% foa_3d                []          plot cones in 3D.
%                                       figure or axes. if 0, make new fig. if empty, don't plot.
% foa_pie               []          plot pie chart of how cones were classified.
% foa_2d                []          plot cones in projection space, as specified in cone_remap.
% foa_hist              []          histogram of sensitivities of L and M cones along the line spearating them.
% cone_remap           	[]          struct specifying a map to the default 2D color space: R/(R+G) vs B/(B+G)
%                                   structure with fields:
%                                       fcn       - function mapping from 3D to 2D to plot RGB from all cones
%                                       x_caption - caption for x axis
%                                       y_caption - caption for y axis
% dot_size              5           size of dots in plot
%
%
%
% 2008-10 gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('foa_3d', []);
p.addParamValue('foa_2d', []);
p.addParamValue('foa_pie', []);
p.addParamValue('foa_hist', []);
p.addParamValue('cone_remap', ...
    struct('fcn',@(x)([x(:,1)./x(:,2) x(:,3)./x(:,2)]),...
    'x_caption','red/green','y_caption','blue/green'));
p.addParamValue('dot_size', 5);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;



% get variables from datarun
cone_rgb = datarun.cones.rgb;
rgb_expected = cone_rgb_expected(datarun);
cone_index.L = find(datarun.cones.types == 'L');
cone_index.M = find(datarun.cones.types == 'M');
cone_index.S = find(datarun.cones.types == 'S');
cone_index.U = find(datarun.cones.types == 'U');
num_cones = size(datarun.cones.rgb,1);



% set up plot axes
pa_3d = set_up_fig_or_axes(params.foa_3d);
pa_2d = set_up_fig_or_axes(params.foa_2d);
pa_pie = set_up_fig_or_axes(params.foa_pie);
pa_hist = set_up_fig_or_axes(params.foa_hist);


% plot cone colors in 2D remapped space
if ~isempty(pa_2d) && ~isempty(params.cone_remap) && isfield(params.cone_remap,'fcn')
    
    % select axes
    axes(pa_2d)

    % map 3D cone weights to 2D
    remapped_rgb = params.cone_remap.fcn(cone_rgb(1:num_cones,1:3));
    
    % plot each cone type in a different color
    plot(remapped_rgb(cone_index.L,1),remapped_rgb(cone_index.L,2),'.r','MarkerSize',params.dot_size);hold on
    plot(remapped_rgb(cone_index.M,1),remapped_rgb(cone_index.M,2),'.g','MarkerSize',params.dot_size)
    plot(remapped_rgb(cone_index.S,1),remapped_rgb(cone_index.S,2),'.b','MarkerSize',params.dot_size)
    plot(remapped_rgb(cone_index.U,1),remapped_rgb(cone_index.U,2),'.k','MarkerSize',params.dot_size)
    
    if isfield(params.cone_remap,'x_caption') && isfield(params.cone_remap,'y_caption')
        xlabel(params.cone_remap.x_caption)
        ylabel(params.cone_remap.y_caption)
    end

    % plot expected gun strengths

    % get L and M center in remapped coordinates
    L_ctr = params.cone_remap.fcn(rgb_expected.L);
    M_ctr = params.cone_remap.fcn(rgb_expected.M);
    % add them to the plot
    hold on
    plot(L_ctr(1),L_ctr(2),'+k','MarkerSize',25)
    plot(M_ctr(1),M_ctr(2),'+k','MarkerSize',25)

    % set axes
    set(gca,'XLim',[-.5 1.5],'YLim',[-.5 1.5])

    drawnow
end


% plot cone RGB values in 3D
if ~isempty(pa_3d)
    
    % select axes
    axes(pa_3d)
    
    % enter rgb triplets into easier variables
    rd = cone_rgb(:,1);
    gn = cone_rgb(:,2);
    be = cone_rgb(:,3);

    % put cone weights into easier variables
    Ln = rgb_expected.L/norm(rgb_expected.L);
    Mn = rgb_expected.M/norm(rgb_expected.M);
    Sn = rgb_expected.S/norm(rgb_expected.S);

    % plot points
    plot3(rd(cone_index.L),gn(cone_index.L),be(cone_index.L),'.r','MarkerSize',params.dot_size) % L cones
    hold on
    plot3(rd(cone_index.M),gn(cone_index.M),be(cone_index.M),'.g','MarkerSize',params.dot_size) % M cones
    plot3(rd(cone_index.S),gn(cone_index.S),be(cone_index.S),'.b','MarkerSize',params.dot_size) % S cones
    plot3(rd(cone_index.U),gn(cone_index.U),be(cone_index.U),'.k','MarkerSize',params.dot_size) % unsure cones
    
    % plot ideal cone lines
    scl = max([get(gca,'XLim') get(gca,'YLim') get(gca,'ZLim')]);
    plot3([0 Ln(1)]*scl,[0 Ln(2)]*scl,[0 Ln(3)]*scl,'r',...  red ideal
        [0 Mn(1)]*scl,[0 Mn(2)]*scl,[0 Mn(3)]*scl,'g',...  green ideal
        [0 Sn(1)]*scl,[0 Sn(2)]*scl,[0 Sn(3)]*scl,'b')     %blue ideal
    
    % label axes
    xlabel('red');ylabel('green');zlabel('blue');
    
    % choose a better view
    set(gca,'CameraPosition',[478 1870 2370])

end



% plot pie chart showing number of each type
if ~isempty(pa_pie)
    
    % select axes
    axes(pa_pie)
    
    pie([length(cone_index.L) length(cone_index.U) length(cone_index.M) length(cone_index.S)],...
        {sprintf('L\n%d',length(cone_index.L)),sprintf('L or M\n%d',length(cone_index.U)),...
        sprintf('M\n%d',length(cone_index.M)),sprintf('S\n%d',length(cone_index.S))})
    colormap([1 0 0; .5 .5 .5; 0 1 0; 0 0 1])
end



% plot histogram showing separation
if ~isempty(pa_hist)
    cone_rgb_renorm = normalize_to_unit_sphere(cone_rgb);
    

    % get centroids of L and M cones
    % these will determine the projection axis
    % only if there is at least one cone of the type
    if ~isempty(cone_index.L)
        ctr_L_rgb = mean(cone_rgb_renorm(cone_index.L,:));
    else
        ctr_L_rgb = [1 1 1];
    end
    if ~isempty(cone_index.M)
        ctr_M_rgb = mean(cone_rgb_renorm(cone_index.M,:));
    else
        ctr_M_rgb = [.5 .5 .5];
    end

    % get overall centroid
    ctr_LM_rgb = (ctr_L_rgb + ctr_M_rgb)/2;

    % re-center L, M centroids so that overall centroid is at the origin
    ctr_L = ctr_L_rgb - ctr_LM_rgb;
    ctr_M = ctr_M_rgb - ctr_LM_rgb;

    % compute projection axis
    proj_axis = ctr_L - ctr_M;
    
    % function to project rgb values
    project = @(x)(x-repmat(ctr_LM_rgb,size(x,1),1))*proj_axis';
    
    % project cone rgb values
    %cone_proj = (cone_rgb_renorm - repmat(ctr_LM,num_cones,1)) * proj_axis';
    cone_proj = project(cone_rgb_renorm);
    
    % separate L and M
    Lp = cone_proj(cone_index.L,:);
    Mp = cone_proj(cone_index.M,:);
    
    % set good bins
    proj_radius = 4*(project(ctr_L_rgb) - project(ctr_M_rgb));
    bins = proj_radius*[-31:31]/62;
    
    % get histogram
    % L cones
    [nL,x] = hist(cone_proj(cone_index.L),bins);
    % M cones
    [nM,x] = hist(cone_proj(cone_index.M),bins);
    % unknown cones
    [nU,x] = hist(cone_proj(cone_index.U),bins);
    
    % plot it
    axes(pa_hist); hold on
    %hL = bar(x(2:end-1),nL(2:end-1),'r')
    %hM = bar(x(2:end-1),nM(2:end-1),'g')
    h = bar(x(2:end-1),[nU(2:end-1); nL(2:end-1); nM(2:end-1)]','stacked');
    set(h(1),'FaceColor',[.5 .5 .5],'LineStyle','none')
    set(h(2),'FaceColor','r','LineStyle','none')
    set(h(3),'FaceColor','g','LineStyle','none')
    xlim([x(2) x(end-1)])
    
    % add ideal cone lines
    plot([1 1]*project(ctr_L_rgb),ylim,'-','LineWidth',1,'Color',[.5 0 0])
    plot([1 1]*project(ctr_M_rgb),ylim,'-','LineWidth',1,'Color',[0 .5 0])
    
    
end


