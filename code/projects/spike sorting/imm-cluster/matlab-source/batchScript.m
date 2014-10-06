range = 1:256;
itName = 'test1';
prjPath = '/blah/blah.prj';
oPath = 'blah;';
noOutput = true;


params = config(range, itName, prjPath, oPath, noOutput);
ImmCluster(params)

exit;
