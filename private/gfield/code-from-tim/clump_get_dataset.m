% clump_get_dataset
%
% get dataset paths. examples:
% call clump_getDataset(0)             to get all the locations
% call clump_getDataset(1:4)           to get datasets 1 through 4
% call clump_getDataset([],'plantain') to get plantain

function [cd, name] = clump_get_dataset(rng, nameString)

% define locations for text files

name{18} = 'simulation';
cd{18} = '/marte/snle/home/gfield/Desktop/cone-simulations/scale-5-clumping/cones.txt';

name{1}  = 'blueberry';
cd{1} = '/marte/snle/lab/Experiments/Array/Shared/one/2008-08-26-2_data001_data001-fit 0.75/cones.txt';
name{2}  = 'grapes';
cd{2} = '/marte/snle/lab/Experiments/Array/Shared/one/2007-03-27-2_data014_data014_data014-fit 0.75/cones.txt';
name{3}  = 'kiwi';
cd{3} = '/marte/snle/lab/Experiments/Array/Shared/one/2008-05-13-3_data006_data006-fit 0.75/cones.txt';
name{4}  = 'mango';
cd{4} = '/marte/snle/lab/Experiments/Array/Shared/one/2008-04-30-2_data004_data004_data004-fit 0.75/cones.txt';
name{5}  = 'plantain';
cd{5} = '/marte/snle/lab/Experiments/Array/Shared/one/2008-08-27-5_data003_data003_data003-fit 0.75/cones.txt';
name{6}  = 'plum';
cd{6} = '/marte/snle/lab/Experiments/Array/Shared/one/2008-04-22-5_data006_data006-fit 0.75/cones.txt';
name{7}  = 'pomegranate';
cd{7} = '/marte/snle/lab/Experiments/Array/Shared/one/2007-08-21-1_data003_data003-fit 0.75/cones.txt';
name{8}  = 'raspberry';
cd{8} = '/marte/snle/lab/Experiments/Array/Shared/one/2006-06-12-0_data003_data003-fit 0.75/cones.txt';
name{8}  = 'butterfly';
cd{9} = '/marte/snle/lab/Experiments/Array/Shared/one/2008-12-12-1_data005_data005-fit 0.75/cones.txt';
name{10}  = 'peach';
cd{10} = '/marte/snle/lab/Experiments/Array/Shared/one/2008-08-27-0_data001-s5661-s9260_data001-s5661-s9260-fit 0.75/cones.txt';
name{11} = 'simulation-15';
cd{11} = '/marte/snle/lab/Experiments/Array/Shared/one/simulations/simulation-balanced-15clumped/cones.txt';
name{12} = 'simulation-30';
cd{12} = '/marte/snle/lab/Experiments/Array/Shared/one/simulations/simulation-balanced-30clumped/cones.txt';
name{13} = 'simulation-90';
cd{13} = '/marte/snle/lab/Experiments/Array/Shared/one/simulations/simulation-balanced-90clumped/cones.txt';
name{14} = 'peach-manual';
cd{14} = '/marte/snle/lab/Experiments/Array/Shared/one/2008-08-27-0_data001-s5661-s9260_data001-s5661-s9260-manual/cones.txt';
name{15} = 'plum-manual';
cd{15} = '/marte/snle/lab/Experiments/Array/Shared/one/2008-04-22-5_data006_data006-manual/cones.txt';
name{16} = 'plantain-manual';
cd{16} = '/marte/snle/lab/Experiments/Array/Shared/one/2008-08-27-5_data003_data003_data003-manual/cones.txt';
name{17} = 'pomegranate-manual';
cd{17} = '/marte/snle/lab/Experiments/Array/Shared/one/2007-08-21-1_data003_data003-manual2/cones.txt';

% return the values we want to use
if nargin < 2
    if rng ~= 0
        cd = {cd{rng}}; name = {name{rng}};
    end
else
    for ii = 1:length(name)
        if isequal(name{ii},nameString)
            cd{1} = cd{ii}; name{1} = name{ii};
            return;
        end
    end
end