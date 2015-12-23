function [LMS_paths, LMS_names, cone_paths] = get_LMS_paths(quality_threshold, varargin)
% Returns the paths and data set names for dataruns used in LMS manuscript

% usage: [LMS_paths, LMS_names] = get_LMS_paths(quality_threshold)

% required arguments:
%          quality_threshold        must be 'high', 'medium', 'low'
%                                   low returns paths to all data 
%                                   medium returns just the paths to good data and OK data
%                                   high returns just the paths to good data
%
p = inputParser;

% specify list of optional parameters
p.addParamValue('cone_finding', 'standard');
p.parse(varargin{:});

cone_finding = p.Results.cone_finding;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% construct database of data set names, paths, and assigned qualities
% blueberry
lms(1).name = 'blueberry';
lms(1).folder = '/snle/lab/Experiments/Array/Analysis/2008-08-26-2/data001-s6369-s9551/data001-s6369-s9551';
lms(1).quality = 2;

% peach
lms(2).name = 'peach';
lms(2).folder = '/snle/lab/Experiments/Array/Analysis/2008-08-27-0/data001-0s-2400s/data001/data001';
lms(2).quality = 2;

% plantain
lms(3).name = 'plantain';
lms(3).folder = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
lms(3).quality = 2;

% kiwi
lms(4).name = 'kiwi';
lms(4).folder = '/snle/lab/Experiments/Array/Analysis/2008-05-13-3/data006/data006';
lms(4).quality = 2;

% grapes
lms(5).name = 'grapes';
lms(5).folder = '/snle/lab/Experiments/Array/Analysis/2007-03-27-2/data014/data014/data014';
lms(5).quality = 2;

% apricot
lms(6).name = 'apricot';
lms(6).folder = '/snle/lab/Experiments/Array/Analysis/2009-04-13-5/data005-s3600-s7200/data005/data005';
lms(6).quality = 2;

% apple
lms(7).name = 'apple';
lms(7).folder = '/snle/lab/Experiments/Array/Analysis/2010-03-05-2/data013/data013';
lms(7).quality = 2;

% plantain % cones estimated from peek frame, no projection along time course
% lms(8).name = 'plantain-2';
% lms(8).folder = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
% lms(8).quality = 2;


lms(8).name = 'butterfly';
lms(8).folder = '/snle/lab/Experiments/Array/Analysis/2008-12-12-1/data005/data005/data005';
lms(8).cone_path = 'butterfly';
lms(8).quality = 1;

switch cone_finding
    case 'standard';  
        lms(5).cone_path = '2007-03-27-2_data014_data014_data014-bayes-msf_25.00--standard';
        lms(4).cone_path = '2008-05-13-3_data006_data006-bayes-msf_85.00--standard';
        lms(1).cone_path = '2008-08-26-2_data001-s6369-s9551_data001-s6369-s9551-bayes-msf_20.00--standard';
        lms(2).cone_path = '2008-08-27-0_data001-0s-2400s_data001_data001-bayes-msf_20.00--standard';
        lms(3).cone_path = '2008-08-27-5_data003_data003_data003-bayes-msf_70.00--standard';
        lms(6).cone_path = '2009-04-13-5_data005-s3600-s7200_data005_data005-bayes-msf_20.00--standard';
        lms(7).cone_path = '2010-03-05-2_rf-13-apple-gf-bayes-msf_25.00--standard';
        %lms(8).cone_path = '2008-08-27-5_data003_data003_data003-bayes-msf_70.00-70-test';
        lms(8).cone_path = '2008-12-12-1_data005_data005-fit_0.75';

    
    case 'conservative';  
        lms(5).cone_path = '2007-03-27-2_data014_data014-bayes-msf_40.00--conservative';
        lms(4).cone_path = '2008-05-13-3_data006_data006-bayes-msf_140.00--conservative';    
        lms(1).cone_path = '2008-08-26-2_data001-s6369-s9551_data001-s6369-s9551-bayes-msf_40.00--conservative';
        lms(2).cone_path = '2008-08-27-0_data001-0s-2400s_data001_data001-bayes-msf_40.00--conservative';
        lms(3).cone_path = '2008-08-27-5_data003_data003_data003-bayes-msf_160.00--conservative';
        lms(6).cone_path = '2009-04-13-5_data005-s3600-s7200_data005_data005-bayes-msf_40.00--conservative';

    case 'liberal';
        lms(5).cone_path = '2007-03-27-2_data014_data014-bayes-msf_10.00--liberal';
        lms(4).cone_path = '2008-05-13-3_data006_data006-bayes-msf_70.00--liberal';
        lms(1).cone_path = '2008-08-26-2_data001-s6369-s9551_data001-s6369-s9551-bayes-msf_10.00--liberal';
        lms(2).cone_path = '2008-08-27-0_data001-0s-2400s_data001_data001-bayes-msf_10.00--liberal';
        lms(3).cone_path = '2008-08-27-5_data003_data003_data003-bayes-msf_40.00--liberal';
        lms(6).cone_path = '2009-04-13-5_data005-s3600-s7200_data005_data005-bayes-msf_10.00--liberal';
    
    case 'no_on_midget'; 
        lms(5).cone_path = '2007-03-27-2_data014_data014-bayes-msf_25.00--no_ON_midget';
        lms(4).cone_path = '2008-05-13-3_data006_data006-bayes-msf_85.00--no_ON_midget';
        lms(1).cone_path = '2008-08-26-2_data001-s6369-s9551_data001-s6369-s9551-bayes-msf_10.00--no_ON_midget';
        lms(2).cone_path = '2008-08-27-0_data001-0s-2400s_data001_data001-bayes-msf_15.00--no_ON_midget';
        lms(3).cone_path = '2008-08-27-5_data003_data003_data003-bayes-msf_85.00--no_ON_midget';
        lms(6).cone_path = '2009-04-13-5_data005-s3600-s7200_data005_data005-bayes-msf_15.00--no_ON_midget';

    case 'no_off_midget';
        lms(5).cone_path = '2007-03-27-2_data014_data014-bayes-msf_15.00--no_OFF_midget';
        lms(4).cone_path = '2008-05-13-3_data006_data006-bayes-msf_70.00--no_OFF_midget';
        lms(1).cone_path = '2008-08-26-2_data001-s6369-s9551_data001-s6369-s9551-bayes-msf_10.00--no_OFF_midget';
        lms(2).cone_path = '2008-08-27-0_data001-0s-2400s_data001_data001-bayes-msf_25.00--no_OFF_midget';
        lms(3).cone_path = '2008-08-27-5_data003_data003_data003-bayes-msf_85.00--no_OFF_midget';
        lms(6).cone_path = '2009-04-13-5_data005-s3600-s7200_data005_data005-bayes-msf_10.00--no_OFF_midget';
end

%  
% lms(8).name = 'plum';
% lms(8).folder = '/snle/lab/Experiments/Array/Analysis/2008-04-22-5/data006/data006';
% lms(8).cone_path = 'plum';
% lms(8).quality = 1;
% 
% lms(9).name = 'cherry';
% lms(9).folder = '/snle/lab/Experiments/Array/Analysis/2009-02-28-0/data006/data006-0s-3600s/data006/data006';
% lms(9).cone_path = 'cherry';
% lms(9).quality = 1;
% 
% lms(10).name = 'mango';
% lms(10).folder = '/snle/lab/Experiments/Array/Analysis/2008-04-30-2/data004/data004/data004';
% lms(10).cone_path = 'mango';
% lms(10).quality = 1;
% 
% lms(11).name = 'pomegranate';
% lms(11).folder = '/snle/lab/Experiments/Array/Analysis/2007-08-21-1/data003/data003';
% lms(11).cone_path = 'pomegranate';
% lms(11).quality = 0;
% 
% lms(12).name = 'cherimoya';
% lms(12).folder = '/snle/lab/Experiments/Array/Analysis/2008-03-25-3/data002/data002/data002';
% lms(12).cone_path = 'cherimoya';
% lms(12).quality = 0;
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_datasets = length(lms);

if strcmp(quality_threshold, 'high');
    for dset = 1:num_datasets
        if lms(dset).quality > 1
            dset_keeper(dset) = 1;
        else
            dset_keeper(dset) = 0;
        end
    end
elseif strcmp(quality_threshold, 'medium');
    for dset = 1:num_datasets
        if lms(dset).quality > 0 
            dset_keeper(dset) = 1;
        else
            dset_keeper(dset) = 0;
        end
    end
elseif strcmp(quality_threshold, 'low');
    dset_keeper(1:num_datasets) = 1;
else
    error('Incorrect argument: quality threshold must be high medium, or low')
end


% make new struct that keeps only specified data sets
keeper_indices = find(dset_keeper == 1);
kept_lms = lms(keeper_indices);

%initalize outputs
LMS_paths = cell(1,length(keeper_indices));
LMS_names = cell(1,length(keeper_indices));
cone_paths = cell(1,length(keeper_indices));

[LMS_paths{:}] = deal(kept_lms.folder);
[LMS_names{:}] = deal(kept_lms.name);
[cone_paths{:}] = deal(kept_lms.cone_path);





