function trakem2_area_scatter(project)

hold on;

somas = project.project_tree.ganglion_cell_layer{1}.soma;
soma_nucleus_mismatch_size = 0;
nucleus_bigger_than_soma_area = 0;

for i = 1:length(somas)
    soma = somas{i};
    soma_max_area = get_max_area(soma.area_list);
    
    nucleus = soma.nucleus{1};
    nucleus_max_area = get_max_area(nucleus.area_list);
    
    if( length( soma_max_area ) == length( nucleus_max_area ) ) 
        if( nucleus_max_area > soma_max_area )
           nucleus_bigger_than_soma_area =  nucleus_bigger_than_soma_area + 1;
        end
        
        plot(soma_max_area, nucleus_max_area, 'o');
        
    else
        soma_nucleus_mismatch_size = soma_nucleus_mismatch_size + 1;
    end
 
end
disp('soma nucleus mismatch number: ');
disp(soma_nucleus_mismatch_size);

disp( 'number of nuclei bigger than soma: ' );
disp( nucleus_bigger_than_soma_area);

function max_area = get_max_area(area_lists)
areas = [];
for j = 1:length(area_lists)
    al = area_lists{j};
    if isfield (al, 'layers' )
        layers = al.layers;
        for k = 1:length(layers)
            layer = layers(k);
            areas = [areas layer.area];
        end
    end
end
max_area = max(areas);