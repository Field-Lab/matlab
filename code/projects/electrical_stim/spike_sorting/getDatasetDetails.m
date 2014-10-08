function details = getDatasetDetails(dataSetID)

%% representative datasets

switch(dataSetID)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1
        details.savePath = '/snle/home/lhruby/Documents/projects/data_analysis/automated/analysis/rep_datasets/2008-08-26-0-data006_pw50/';
        %details.dataPath = '/Volumes/Lee/Analysis/Lauren/2008-08-26-0/data006';
        %details.artifactPath = '/Volumes/Lee/Analysis/Lauren/2008-08-26-0/data011';
        details.dataPath = '/snle/lab/Experiments/Array/Analysis/2008-08-26-0/data006';
        details.artifactPath = '/snle/lab/Experiments/Array/Analysis/2008-08-26-0/data011';

        details.neuronID = 151;
        details.neuronIDsFull = [2 182 602 183 31 106 108 151 109];
        
        %details.pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-08-26-0/data005/data005.ei';
        details.pathToEi = '/snle/lab/Experiments/Array/Analysis/2008-08-26-0/data005/data005.ei';
        
        details.centerChannel = 7;
        details.channelRadius = 1;

        details.patternNumbers = 7;
        details.movieNumbers = 1:3:76;
        details.artifactPatternNumbers = 7;
        details.artifactMovieNumbers = 1:3:76;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2
        details.savePath  = '/snle/home/lhruby/Documents/projects/data_analysis/automated/analysis/rep_datasets/2008-08-26-0-data006_pw100/';
        %details.dataPath = '/Volumes/Lee/Analysis/Lauren/2008-08-26-0/data006';
        %details.artifactPath = '/Volumes/Lee/Analysis/Lauren/2008-08-26-0/data011';
        details.dataPath = '/snle/lab/Experiments/Array/Analysis/2008-08-26-0/data006';
        details.artifactPath = '/snle/lab/Experiments/Array/Analysis/2008-08-26-0/data011';
        
        details.neuronID = 151;
        details.neuronIDsFull = [2 182 602 183 31 106 108 151 109];
        %details.pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-08-26-0/data005/data005.ei';
        details.pathToEi = '/snle/lab/Experiments/Array/Analysis/2008-08-26-0/data005/data005.ei';
        
        details.centerChannel = 7;
        details.channelRadius = 1;

        details.patternNumbers = 7;
        details.movieNumbers = 2:3:74;
        details.artifactPatternNumbers = 7;
        details.artifactMovieNumbers = 2:3:74;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 3
        details.savePath = '/snle/home/lhruby/Documents/projects/data_analysis/automated/analysis/rep_datasets/2008-08-27-2-data002/';
        %details.dataPath = '/Volumes/Lee/Analysis/Lauren/2008-08-27-2/data002';
        details.artifactPath = '';
        details.dataPath = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data002';

        details.neuronID = 91;
        details.neuronIDFull = [91 139 108 902 18 19];
        %details.pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-08-27-2/data001-lh/data001-lh.ei';
        details.pathToEi = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data001-lh/data001-lh.ei';
        
        details.centerChannel = 3;
        details.channelRadius = 1;

        details.patternNumbers = 3;
        details.movieNumbers = 1:26;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 4
        details.savePath = '/snle/home/lhruby/Documents/projects/data_analysis/automated/analysis/rep_datasets/2008-08-27-4-data001_p3/';
        %details.dataPath = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data001';
        %details.artifactPath = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data013';
        details.dataPath = '/snle/lab/Experiments/Array/Analysis/2008-08-27-4/data001';
        details.artifactPath = '/snle/lab/Experiments/Array/Analysis/2008-08-27-4/data013';

        details.neuronID = [887 32];
        details.neuronIDFull = [887 32 79 110];
        %details.pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data000-lh/data000-lh.ei';
        details.pathToEi = '/snle/lab/Experiments/Array/Analysis/2008-08-27-4/data000-lh/data000-lh.ei';
        
        details.centerChannel = 3;
        details.channelRadius = 1;

        details.patternNumbers = 3;
        details.movieNumbers = 1:25;
        details.artifactPatternNumbers = 3;
        details.artifactMovieNumbers = 2:3:74;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 5
        details.savePath = '/snle/home/lhruby/Documents/projects/data_analysis/automated/analysis/rep_datasets/2008-08-27-4-data001_p8/';
        %details.dataPath = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data001';
        %details.artifactPath = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data013';
        details.dataPath = '/snle/lab/Experiments/Array/Analysis/2008-08-27-4/data001';
        details.artifactPath = '/snle/lab/Experiments/Array/Analysis/2008-08-27-4/data001';
        
        
        
        details.neuronID = 107;
        details.neuronIDFull = [107 106 108 137];
        %details.pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data000-lh/data000-lh.ei';
        details.pathToEi = '/snle/lab/Experiments/Array/Analysis/2008-08-27-4/data000-lh/data000-lh.ei';
        
        details.centerChannel = 8;
        details.channelRadius = 1;

        details.patternNumbers = 8;
        details.movieNumbers = 1:25;
        details.artifactPatternNumbers = 8;
        details.artifactMovieNumbers = 2:3:74;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 6
        details.savePath = '/netapp/snle/home/lhruby/Documents/projects/data_analysis/automated/analysis/rep_datasets/2008-11-10-2-data001/';
        %details.dataPath = '/Volumes/Lee/Analysis/Lauren/2008-11-10-2/data001';
        details.artifactPath = '';
        details.dataPath = '/snle/lab/Experiments/Array/Analysis/2008-11-10-2/data001';
        
        details.neuronID = [92 106 228];
        details.neuronIDFull = [92 106 228 156];
        %details.pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-11-10-2/data000/data000.ei';
        details.pathToEi = '/snle/lab/Experiments/Array/Analysis/2008-11-10-2/data000/data000.ei';

        details.centerChannel = 8;
        details.channelRadius = 1;

        details.patternNumbers = 8;
        details.movieNumbers = 1:26;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% has TTX artifact recordings, but they are incomplete
    case 7
        details.savePath = '/snle/home/lhruby/Documents/projects/data_analysis/automated/analysis/rep_datasets/2008-11-10-3-data003_pw50/';
        %details.dataPath = '/Volumes/Lee/Analysis/Lauren/2008-11-10-3/data003';
        %details.artifactPath = '/Volumes/Lee/Analysis/Lauren/2008-11-10-3/data020'; %%higher amp artifacts in data021
        details.artifactPath = '';
        details.dataPath = '/snle/lab/Experiments/Array/Analysis/2008-11-10-3/data003';
        
        details.neuronID = [648 677];
        details.neuronIDFull = [648 677 661 511 544 634 767 769];
        %details.pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-11-10-3/data004/data004.ei';
        details.pathToEi = '/snle/lab/Experiments/Array/Analysis/2008-11-10-3/data004/data004.ei';

        details.centerChannel = 46;
        details.channelRadius = 1;

        details.patternNumbers = 46;
        details.movieNumbers = 1:3:97;
        %details.artifactPatternNumbers = 46;
        %details.artifactMovieNumbers1 = 1:3:;
        %details.artifactMovieNumbers2 = :3:97;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% has TTX artifact recordings, but they are incomplete
    case 8
        details.savePath = '/snle/home/lhruby/Documents/projects/data_analysis/automated/analysis/rep_datasets/2008-11-10-3-data003_pw100/';
        %details.dataPath = '/Volumes/Lee/Analysis/Lauren/2008-11-10-3/data003';
        %details.artifactPath = '/Volumes/Lee/Analysis/Lauren/2008-11-10-3/data020'; %%higher amp artifacts in data021
        details.artifactPath = '';
        details.dataPath = '/snle/lab/Experiments/Array/Analysis/2008-11-10-3/data003';
        
        details.neuronID = [648 677];
        details.neuronIDFull = [648 677 661 511 544 634 767 769];
        %details.pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-11-10-3/data004/data004.ei';
        details.pathToEi = '/snle/lab/Experiments/Array/Analysis/2008-11-10-3/data004/data004.ei';
        
        details.centerChannel = 46;
        details.channelRadius = 1;

        details.patternNumbers = 46;
        details.movieNumbers = 2:3:95;
        %details.artifactPatternNumbers = 46;
        %details.artifactMovieNumbers1 = 2:3:;
        %details.artifactMovieNumbers2 = :3:95;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 9
        details.savePath = '/snle/home/lhruby/Documents/projects/data_analysis/automated/analysis/rep_datasets/2008-11-12-3-data001/';
        %details.dataPath = '/Volumes/Lee/Analysis/Lauren/2008-11-12-3/data001';
        %details.artifactPath = '/Volumes/Lee/Analysis/Lauren/2008-11-12-3/data005';
        details.dataPath = '/snle/lab/Experiments/Array/Analysis/2008-11-12-3/data001';
        details.artifactPath = '/snle/lab/Experiments/Array/Analysis/2008-11-12-3/data005';
        
        details.neuronID = [904 931];
        details.neuronIDFull = [873 904 931];
        %details.pathToEi = '/Volumes/Lee/Analysis/Lauren/2008-11-12-3/data000/data000.ei';
        details.pathToEi = '/snle/lab/Experiments/Array/Analysis/2008-11-12-3/data000/data000.ei';
        
        details.centerChannel = 62;
        details.channelRadius = 1;

        details.patternNumbers = 62;
        details.movieNumbers = 1:26;
        details.artifactPatternNumbers = 62;
        details.artifactMovieNumbers = 1:26;
        
    otherwise
        error('No dataset match.')
end