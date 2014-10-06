function project_tree = parse_jproject_tree(jpt)
root = jpt.getRoot();
project_tree = parse_jproject_tree_node(root);


function parsed = parse_jproject_tree_node(node)
parsed = struct();
parsed.jnode = node;
UO = node.getUserObject();
parsed.title = char(UO.getTitle());
parsed.type  = char(UO.getType());

numchildren = node.getChildCount();
parsed.children = cell(numchildren, 1);
for i = 1:numchildren
    child = node.getChildAt(i-1);
    parsed_child = parse_jproject_tree_node(child);
    parsed.children{i} = parsed_child;
    
    childtype = parsed_child.type;
    if isfield(parsed, childtype)
        parsed.(childtype){end+1} = parsed_child;
    else
        parsed.(childtype) = {parsed_child};
    end
end

if strcmp(parsed.type, 'area_list')
    al = UO.getObject();
    area = al.getArea();
    bounds = area.getBounds();
    parsed.bounds = [bounds.getX bounds.getY];
    parsed.width  = bounds.getWidth;
    parsed.height = bounds.getHeight;
    
    jlayerids = al.getLayerIds;
    for i = 1:jlayerids.size()
        layerid = jlayerids.get(i-1);
        parsed.layers(i).id = layerid;

        jpaths = al.getPaths(layerid);
        parsed.layers(i).paths =  parse_jpaths(jpaths);
        
        parsed.layers(i).masks = paths2masks(parsed.layers(i).paths, parsed.width, parsed.height);
        parsed.layers(i).area = total_area(parsed.layers(i).paths);
    end    
end


function paths = parse_jpaths(jpaths)
paths = cell(jpaths.size(), 1);
for i = 1:length(paths)
    jpath = jpaths.get(i-1);
    paths{i} = parse_jpath(jpath);
end


function path = parse_jpath(jpath)
path = zeros(jpath.size(), 2);
for i = 1:jpath.size()
    point = jpath.get(i-1);
    path(i,1) = point.getX;
    path(i,2) = point.getY;
end


function masks = paths2masks(paths, width, height)
masks = cell(length(paths), 1);
for i = 1:length(masks)
    path = paths{i};
    if ~isempty(path)
        masks{i} = poly2mask(path(:,1), path(:,2), height, width)';
    end
end


function area = total_area(paths)
area = 0;
for i = 1:length(paths)
    area = area + polyarea(paths{i}(:,1), paths{i}(:,2));
end