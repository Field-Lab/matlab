function dsInfo =  dataset_list_axon_directions_StanfordDirs()

%list of ei files with corresponding neuron ids that have decent axon
%signals for tracing

x = 0;

x = x+1;
dsInfo(x).eiPath = '/Volumes/Analysis/2008-08-26-0/data007-NW/data007-NW.ei';
dsInfo(x).axonTraceCells = [334 406 511 526 541 557 560 559 588];

x = x+1;
dsInfo(x).eiPath = '/Volumes/Analysis/2008-08-27-4/data004-NW/data004-NW.ei';
dsInfo(x).axonTraceCells = [18 182 31 106 151 871 61 183 196 933];

x = x+1;
dsInfo(x).eiPath = '/Volumes/Analysis/2008-11-12-3/data000/data000.ei';
dsInfo(x).axonTraceCells = [904 812 931];

x = x+1;
dsInfo(x).eiPath = '/Volumes/Analysis/2011-01-11-0/data030/data030.ei';
dsInfo(x).axonTraceCells = [34 871 33 62 799 814 877 1 932 767 874];

x = x+1;
dsInfo(x).eiPath = '/Volumes/Analysis/2011-05-11-8/data001/data001.ei';
dsInfo(x).axonTraceCells = [631 725 811 766 379];

x = x+1;
dsInfo(x).eiPath = '/Volumes/Analysis/2011-07-05-5/data000/data000.ei';
dsInfo(x).axonTraceCells = [316 436 407];

x = x+1;
dsInfo(x).eiPath = '/Volumes/Analysis/2011-07-14-7/data000/data000.ei';
dsInfo(x).axonTraceCells = [406 451 559 586 632 677 781 484 797 143 711 528];

x = x+1;
dsInfo(x).eiPath = '/Volumes/Analysis/2011-08-01-0/data000/data000.ei';
dsInfo(x).axonTraceCells = [1 91 166 181 196 16 857 79 136];

x = x+1;
dsInfo(x).eiPath = '/Volumes/Analysis/2011-08-04-5/data000/data000.ei';
dsInfo(x).axonTraceCells = [766 392 316 728 741 772 949 287 707 736 811 948];

x = x+1;
dsInfo(x).eiPath = '/Volumes/Analysis/2011-12-13-0/data005/data005.ei';
dsInfo(x).axonTraceCells = [544 781 511 631 316 646 691 706];