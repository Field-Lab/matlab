function project = trakem2read(path)
fsloader = ini.trakem2.persistence.FSLoader;
data = fsloader.openFSProject(path, false);
clear fsloader;

% Build project object
jproject = data(2).getObject();
project = struct();
project.jtemplate_tree = ini.trakem2.tree.TemplateTree(jproject, data(1));
project.jproject_tree  = ini.trakem2.tree.ProjectTree( jproject, data(2));
project.jlayer_tree    = ini.trakem2.tree.LayerTree(   jproject, data(3));
clear jproject;

project.project_tree = parse_jproject_tree(project.jproject_tree);
