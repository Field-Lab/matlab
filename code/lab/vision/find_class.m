function jclass = find_class(classname, classloader)
% FIND_CLASS  Get the actual Java class object, rather than an instance
% usage: jclass = find_class(classname, classloader)
%
% I don't know why this is so hard; hopefully there is a better way...
%
% For classes that are easy to instantiate, you can just do: 
%   package.Class().getClass();
%
% But if the constructor is non-trivial or you are trying to get an inner
% class, this may be your only option.  The tricky thing here then is to 
% get a hold of the proper classloader object.  The easiest kludge is to
% pick any trivially constructable object from the same package and use:
%   trivial_object.getClass().getClassLoader();
%
% The default is to do just this with the vision Matlab base class:
%    classloader = edu.ucsc.neurobiology.vision.matlab.Matlab().getClass().getClassLoader();
%
% 2011-09 phli
%

if nargin < 2
    classloader = edu.ucsc.neurobiology.vision.matlab.Matlab().getClass().getClassLoader();
end

jclass = java.lang.Class.forName(classname, true, classloader);