function plot_trakem2(project, min_area, offsets)

if nargin < 2
    min_area = 0;
end

if nargin < 3
    offsets = [0, 0];
end


hold on;

somas = project.project_tree.ganglion_cell_layer{1}.soma;
for i = 1:length(somas)
    soma = somas{i};
    if get_max_area(soma.area_list) < min_area
        continue;
    end
    
    als = soma.area_list;
    for j = 1:length(als)
        al = als{j};
        plot_area_list(al, 'g', offsets);
    end
    
%     nucleus = soma.nucleus{1};
%     als = nucleus.area_list;
%     for j = 1:length(als)
%         plot_area_list(als{j}, 'b');
%     end
end


function plot_area_list(area_list, styles, offsets)
bounds = area_list.bounds;
layers = area_list.layers;
for k = 1:length(layers)
    layer = layers(k);
    paths = layer.paths;
    if ~isempty(paths)
        for l = 1:length(paths)
            path = paths{l};
            if ~isempty(path)
                x = path(:,1) + bounds(1) + offsets(1); % GAH!  Not sure why these are necessary :(
                y = path(:,2) + bounds(2) + offsets(2);
                plot(x,y,styles);
            end
        end
    end
end


function max_area = get_max_area(area_lists)
areas = [];
for j = 1:length(area_lists)
    al = area_lists{j};
    layers = al.layers;
    for k = 1:length(layers)
        layer = layers(k);
        areas = [areas layer.area];
    end
end
max_area = max(areas);