function plot_trakem2_paths(project)

hold on;

somas = project.project_tree.ganglion_cell_layer{1}.soma;
num_soma_plots = 0;

for i = 1:length(somas)
    soma = somas{i};
    als = soma.area_list;
    
    for j = 1:length(als)
        al = als{j};
        plot_area_list(al, 'g');
        num_soma_plots = num_soma_plots + 1;
    end
    
    nucleus = soma.nucleus{1};
    als = nucleus.area_list;
    for j = 1:length(als)
        plot_area_list(als{j}, 'b');
    end
end
disp( num_soma_plots)


function plot_area_list(area_list, styles)
bounds = area_list.bounds;
if isfield( area_list, 'layers')
    layers = area_list.layers;
    for k = 1:length(layers)
        layer = layers(k);
        paths = layer.paths;
        if ~isempty(paths)
            for l = 1:length(paths)
                path = paths{l};
                if ~isempty(path)
                    x = path(:,1) + bounds(1);
                    y = path(:,2) + bounds(2);
                    plot(x,y,styles);
                end
            end
        end
    end
end