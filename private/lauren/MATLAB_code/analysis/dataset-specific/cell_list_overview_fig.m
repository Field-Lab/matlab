function [cellInfo datarun] = cell_list_overview_fig()


%% ON parasol example

cellInfo(1).pathToAnalysis = '/snle/lab/Experiments/Array/Analysis/2011-07-05-0/data002/';
cellInfo(1).id = 527;
cellInfo(1).stimElec = 33;
cellInfo(1).PW = 100;
cellInfo(1).suffix = '_w100';
cellInfo(1).type = 'onPar';
cellInfo(1).gain = 440;

datarun(1).pathToData = '2011-07-05-0/data000/data000';
datarun(1).obvius_fit_path = '/snle/lab/Experiments/Array/Analysis/2011-07-05-0/data000';
datarun(1).cellsToInclude = [78 302 467 527 542 663 798];


%% OFF parasol example

cellInfo(2).pathToAnalysis = '/snle/lab/Experiments/Array/Analysis/2011-10-25-4/data002/';
cellInfo(2).id = 46;
cellInfo(2).stimElec = 1;
cellInfo(2).PW = 100;
cellInfo(2).suffix = '_w100';
cellInfo(2).type = 'offPar';
cellInfo(2).gain = 440;

datarun(2).pathToData = '2011-10-25-4/data000/data000';
datarun(2).obvius_fit_path = '/snle/lab/Experiments/Array/Analysis/2011-10-25-4/data000';
datarun(2).cellsToInclude = [46 198 349 482 647 933];


%% ON midget example

cellInfo(3).pathToAnalysis = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data002/';
cellInfo(3).id = 886;
cellInfo(3).stimElec = 60;
cellInfo(3).PW = 100;
cellInfo(3).suffix = '';
cellInfo(3).type = 'onMidg';
cellInfo(3).gain = 840;

datarun(3).pathToData = ('2008-08-27-2/data001-lh/data001-lh');
datarun(3).obvius_fit_path = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data001-lh';
datarun(3).cellsToInclude = [886 18 857 858 782 872 931];

%% OFF midget example

cellInfo(4).pathToAnalysis = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data029/';
cellInfo(4).id = 680;
cellInfo(4).stimElec = 46;
cellInfo(4).PW = 50;
cellInfo(4).suffix = '_w50';
cellInfo(4).type = 'offMidg';
cellInfo(4).gain = 440;

datarun(4).pathToData = '2011-01-11-0/data030-lh/data030-lh';
datarun(4).obvius_fit_path = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data030-lh';
datarun(4).cellsToInclude = [680 578 624 647 754 799 814];

%% SBC example

cellInfo(5).pathToAnalysis = '/snle/lab/Experiments/Array/Analysis/2010-10-28-2/data003/';
cellInfo(5).id = 18;
cellInfo(5).stimElec = 5;
cellInfo(5).PW = 50;
cellInfo(5).suffix = '_w50';
cellInfo(5).type = 'sbc';
cellInfo(5).gain = 440;

datarun(5).pathToData = '2010-10-28-2/data000/data000';
datarun(5).obvius_fit_path = '/snle/lab/Experiments/Array/Analysis/2010-10-28-2/data000';
datarun(5).cellsToInclude = [18 136 271 481 767];

% extra stuff for SBC only
datarun(5).pathToDataBlur = '2010-10-28-2/data000-blur-0_8/data000';
datarun(5).obvius_fit_path_blur = '/snle/lab/Experiments/Array/Analysis/2010-10-28-2/data000-blur-0_8';
datarun(5).matchedFits = [18 271 481 767]; %cells with qualitatively similar Gaussian fits in blurred and unblurred cases 
%(used to determine STA size correction scale factor)



