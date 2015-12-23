function trakem2_jinit(fijiapppath)

if nargin < 1
    fijiapppath = '/Applications/Fiji.app';
end

javaaddpath([fijiapppath '/jars/ij.jar']);
javaaddpath([fijiapppath '/plugins/TrakEM2_.jar']);
javaaddpath([fijiapppath '/plugins/loci_tools.jar']);
javaaddpath([fijiapppath '/jars/mpicbg.jar']);
javaaddpath([fijiapppath '/plugins/3D_Viewer.jar']);
javaaddpath([fijiapppath '/jars/VectorString.jar']);
javaaddpath([fijiapppath '/jars/imglib.jar']);