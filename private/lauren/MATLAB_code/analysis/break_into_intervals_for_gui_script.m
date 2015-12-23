
clear all

%writePath = '/snle/lab/Experiments/Array/Analysis/2010-11-22-1/data006-filtered/';
writePath = '/Analysis/lauren/2010-11-22-1/data007-filtered/';
%pathToDataFile = '/snle/lab/Experiments/Array/Analysis/2010-11-22-1/data006-filtered/data006-filtered/';
pathToDataFile = '/Analysis/lauren/2010-11-22-1/data007-filtered/data007-filtered/';

interval = 40; %in samples
intsPerMovie = 250;
movieRepsPerMovieOut = 5; %number of original movies repetitions to break up into smaller repetitions of a single output movie number
%movieTimesPath = '/snle/lab/Experiments/Array/Analysis/2010-11-22-1/data006-filtered/data006-movieTimes.mat';
movieTimesPath = '/Analysis/lauren/2010-11-22-1/data007-filtered/data007-movieTimes.mat';
patternNo = 42;
patternStrOrig = '42_f5';
%origAnalysisPath = '/snle/lab/Experiments/Array/Analysis/2010-11-22-1/data006/';
origAnalysisPath = '/Analysis/lauren/2010-11-22-1/data007/';


load(movieTimesPath);

writtenData = break_into_intervals_for_gui_analysis(pathToDataFile, writePath,...
    interval, movieRepsPerMovieOut, movieNumberTimes, patternNo, patternStrOrig, origAnalysisPath);
