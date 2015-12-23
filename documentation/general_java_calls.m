% See COMMON_JAVA_CALLS for Vision specific calls; these are more general
% Matlab/Java interaction patterns

%% Get information on methods including method signatures
methods edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile -full


%% Access classes, inner classes, and enums

ColorType = find_class('edu.ucsc.neurobiology.vision.stimulus.FrameGenerator$ColorType');

% Then we can do nice things like
list_enum(ColorType);
load_enum(ColorType);
INDEPENDENT = get_enum(ColorType, 'INDEPENDENT');

% Actually, if it's an Enum there is a shortcut for just getting values by
% string name:
INDEPENDENT = javaMethod('valueOf', 'edu.ucsc.neurobiology.vision.stimulus.FrameGenerator$ColorType', 'INDEPENDENT');

% Now we can finally instantiate a BinaryFrameGenerator, which required the ColorType!
rng = edu.ucsc.neurobiology.vision.math.RandomJavaV2(seed);
bfg = edu.ucsc.neurobiology.vision.stimulus.BinaryFrameGenerator(100, 100, 0.48, rng, INDEPENDENT);

% (If it's not an Enum and you want to do other things with the class or
% inner class, you can get an instance:)
inner_class.newInstance();


%% Call static inner class method
javaMethod('fromFilename', 'edu.ucsc.neurobiology.vision.io.STACollection$Factory', '/snle/analysis/2011-10-25-9/streaming/data001-0/data001-0.sta')