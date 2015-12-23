function [neuron special_neuron] = cell_list_low_freq_stim()

%note: other datasets with 5 Hz stimulation exist; this lists examples
%of each cell type (4 ON-P, 3 OFF-P, 3 ON-M, 3 OFF-M, 3 SBC)
%
% special_neuron = other OFF-M examples that had to be analyzed differently due to lack of
% spontaneous spikes in long-latency time windows


special_neuron(1).id = 648;
special_neuron(1).stimElec = 44;
special_neuron(1).splitElecRespPath = '/snle/lab/Experiments/Array/Analysis/2010-11-22-1/data006-filtered/';
special_neuron(1).movieTimePath = '/snle/lab/Experiments/Array/Analysis/2010-11-22-1/data006-filtered/data006-movieTimes.mat';
special_neuron(1).type = 'offMidg';
special_neuron(1).originalElecRespPath = '/snle/lab/Experiments/Array/Analysis/2010-11-22-1/data006/elecResp_n648_p44_f5.mat';
special_neuron(1).repStartTimesPath = '/snle/lab/Experiments/Array/Analysis/2010-11-22-1/data006-filtered/repStartTimes.mat';
special_neuron(1).scanElecRespPath = '/snle/lab/Experiments/Array/Analysis/2010-11-22-1/data004/elecResp_n648_p44.mat';
special_neuron(1).PW = 100;

special_neuron(2).id = 621;
special_neuron(2).stimElec = 42;
special_neuron(2).splitElecRespPath = '/snle/lab/Experiments/Array/Analysis/2010-11-22-1/data007-filtered/';
special_neuron(2).movieTimePath = '/snle/lab/Experiments/Array/Analysis/2010-11-22-1/data007-filtered/data007-movieTimes.mat';
special_neuron(2).type = 'offMidg';
special_neuron(2).originalElecRespPath = '/snle/lab/Experiments/Array/Analysis/2010-11-22-1/data007/elecResp_n621_p42_f5.mat';
special_neuron(2).repStartTimesPath = '/snle/lab/Experiments/Array/Analysis/2010-11-22-1/data007-filtered/repStartTimes.mat';
special_neuron(2).scanElecRespPath = '/snle/lab/Experiments/Array/Analysis/2010-11-22-1/data004/elecResp_n621_p42.mat';
special_neuron(2).PW = 100;


z = 0;

z = z+1;
neuron(z).path = '/snle/lab/Experiments/Array/Analysis/2010-10-18-3/data022-filtered/data022-filtered.neurons';
neuron(z).movieTimePath = '/snle/lab/Experiments/Array/Analysis/2010-10-18-3/data022-filtered/data022-movieTimes.mat';
neuron(z).id = 273; %n274 in data004.ei
neuron(z).stimElec = 16;
neuron(z).type = 'onMidg';
neuron(z).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2010-10-18-3/data022/elecResp_n274_p16_f5.mat';
neuron(z).scanElecRespPath = '/snle/lab/Experiments/Array/Analysis/2010-10-18-3/data021/elecResp_n274_p16.mat';
neuron(z).PW = 100;
%for plotting ACF
neuron(z).WN.neuronsPath = '/snle/lab/Experiments/Array/Analysis/2010-10-18-3/data020/data020.neurons';
neuron(z).WN.id = 274;

z = z+1;
neuron(z).path = '/snle/lab/Experiments/Array/Analysis/2010-10-28-2/data006-filtered/data006-filtered.neurons';
neuron(z).movieTimePath = '/snle/lab/Experiments/Array/Analysis/2010-10-28-2/data006-filtered/data006-movieTimes.mat';
neuron(z).id = 64; %18 in data000.ei
neuron(z).stimElec = 5;
neuron(z).type = 'sbc';
neuron(z).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2010-10-28-2/data006/elecResp_n18_p5_f5.mat';
neuron(z).scanElecRespPath = '/snle/lab/Experiments/Array/Analysis/2010-10-28-2/data003/elecResp_n18_p5_w100.mat';
neuron(z).PW = 100;
%for plotting ACF
neuron(z).WN.neuronsPath = '/snle/lab/Experiments/Array/Analysis/2010-10-28-2/data004/data004.neurons';
neuron(z).WN.id = 18;

z = z+1;
neuron(z).path = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data012-filtered/data012-filtered.neurons';
neuron(z).movieTimePath = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data012-filtered/data012-movieTimes.mat';
neuron(z).id = 3; %n947 in data004.ei
neuron(z).stimElec = 58;
neuron(z).type = 'offPar';
neuron(z).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data012/elecResp_n947_p58_f5.mat';
neuron(z).scanElecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data029/elecResp_n888_p58_w100.mat';
neuron(z).PW = 100;
%for plotting ACF
neuron(z).WN.neuronsPath = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data004/data004.neurons';
neuron(z).WN.id = 947;

z = z+1;
neuron(z).path = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data012-filtered/data012-filtered.neurons';
neuron(z).movieTimePath = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data012-filtered/data012-movieTimes.mat';
neuron(z).id = 6; %n2 in data004.ei
neuron(z).stimElec = 58;
neuron(z).type = 'onPar';
neuron(z).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data012/elecResp_n2_p58_f5.mat';
neuron(z).scanElecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data029/elecResp_n2_p58_w100.mat';
neuron(z).PW = 100;
%for plotting ACF
neuron(z).WN.neuronsPath = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data004/data004.neurons';
neuron(z).WN.id = 2;

z = z+1;
neuron(z).path = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data020-filtered/data020-filtered.neurons';
neuron(z).movieTimePath = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data020-filtered/data020-movieTimes.mat';
neuron(z).id = 275; %n230 in data004.ei
neuron(z).type = 'onPar';
neuron(z).stimElec = 20;
neuron(z).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data020/elecResp_n230_p20_f5.mat';
neuron(z).scanElecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data029/elecResp_n257_p20_w100.mat';
neuron(z).PW = 100;
%for plotting ACF
neuron(z).WN.neuronsPath = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data004/data004.neurons';
neuron(z).WN.id = 230;

z = z+1;
neuron(z).path = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data028-filtered/data028-filtered.neurons';
neuron(z).movieTimePath = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data028-filtered/data028-movieTimes.mat';
neuron(z).id = 889; %n875 in data004.ei
neuron(z).stimElec = 61;
neuron(z).type = 'onMidg';
neuron(z).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data028/elecResp_n875_p61_f5.mat';
neuron(z).scanElecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data029/elecResp_n902_p61_w100.mat';
neuron(z).PW = 100;
%for plotting ACF
neuron(z).WN.neuronsPath = '/snle/lab/Experiments/Array/Analysis/2011-01-11-0/data030-lh/data030-lh.neurons';
neuron(z).WN.id = 902;

z = z+1;
neuron(z).path = '/snle/lab/Experiments/Array/Analysis/2011-01-26-3/data004-filtered/data004-filtered.neurons';
neuron(z).movieTimePath = '/snle/lab/Experiments/Array/Analysis/2011-01-26-3/data004-filtered/data004-movieTimes.mat';
neuron(z).id = 785; %n753 in data000.ei
neuron(z).stimElec = 51;
neuron(z).type = 'onPar';
neuron(z).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-01-26-3/data004/elecResp_n753_p51_f5.mat';
neuron(z).scanElecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-01-26-3/data002/elecResp_n753_p51.mat';
neuron(z).PW = 100;
%for plotting ACF
neuron(z).WN.neuronsPath = '/snle/lab/Experiments/Array/Analysis/2011-01-26-3/data000/data000.neurons';
neuron(z).WN.id = 753;

z = z+1;
neuron(z).path = '/snle/lab/Experiments/Array/Analysis/2011-01-26-3/data005-filtered/data005-filtered.neurons';
neuron(z).movieTimePath = '/snle/lab/Experiments/Array/Analysis/2011-01-26-3/data005-filtered/data005-movieTimes.mat';
neuron(z).id = 664; %n616 in data000.ei
neuron(z).stimElec = 45;
neuron(z).type = 'onMidg';
neuron(z).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-01-26-3/data005/elecResp_n616_p45_f5.mat';
neuron(z).scanElecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-01-26-3/data002/elecResp_n616_p45.mat';
neuron(z).PW = 100;
%for plotting ACF
neuron(z).WN.neuronsPath = '/snle/lab/Experiments/Array/Analysis/2011-01-26-3/data000/data000.neurons';
neuron(z).WN.id = 616;

z = z+1;
neuron(z).path = '/snle/lab/Experiments/Array/Analysis/2011-01-26-3/data008-filtered/data008-filtered.neurons';
neuron(z).movieTimePath = '/snle/lab/Experiments/Array/Analysis/2011-01-26-3/data008-filtered/data008-movieTimes.mat';
neuron(z).id = 694; %n482 in data000.ei
neuron(z).stimElec = 36;
neuron(z).type = 'onPar';
neuron(z).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-01-26-3/data008/elecResp_n482_p36_f5.mat';
neuron(z).scanElecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-01-26-3/data002/elecResp_n482_p36.mat';
neuron(z).PW = 100;
%for plotting ACF
neuron(z).WN.neuronsPath = '/snle/lab/Experiments/Array/Analysis/2011-01-26-3/data000/data000.neurons';
neuron(z).WN.id = 482;

% 30-micron array
% z = z+1;
% neuron(z).path = '/snle/lab/Experiments/Array/Analysis/2011-01-26-6/data006-filtered/data006-filtered.neurons';
% neuron(z).movieTimePath = '/snle/lab/Experiments/Array/Analysis/2011-01-26-6/data006-filtered/data006-movieTimes.mat';
% neuron(z).id = 602; %n646 in data000.ei
% neuron(z).stimElec = 36;
% neuron(z).type = 'onMidg';
% neuron(z).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-01-26-6/data006/elecResp_n646_p36_f5.mat';
% neuron(z).scanElecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-01-26-6/data001/elecResp_n646_p36.mat';
% neuron(z).PW = 100;

z = z+1;
neuron(z).path = '/snle/lab/Experiments/Array/Analysis/2011-05-11-2/data011-filtered-cf-lh/data011-filtered-cf-lh.neurons';
neuron(z).movieTimePath = '/snle/lab/Experiments/Array/Analysis/2011-05-11-2/data011-filtered/data011-movieTimes.mat';
neuron(z).id = 592; %n590 in data009.ei
neuron(z).stimElec = 40;
neuron(z).type = 'offPar';
neuron(z).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-05-11-2/data011/elecResp_n590_p40_f5.mat';
neuron(z).scanElecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-05-11-2/data013/elecResp_n590_p40.mat';
neuron(z).PW = 100;
%for plotting ACF
neuron(z).WN.neuronsPath = '/snle/lab/Experiments/Array/Analysis/2011-05-11-2/data009/data009.neurons';
neuron(z).WN.id = 590;

z = z+1;
neuron(z).path = '/snle/lab/Experiments/Array/Analysis/2011-05-11-2/data014-filtered/data014-filtered.neurons';
neuron(z).movieTimePath = '/snle/lab/Experiments/Array/Analysis/2011-05-11-2/data014-filtered/data014-movieTimes.mat';
neuron(z).id = 801; %n797 in data009.ei
neuron(z).stimElec = 47;
neuron(z).type = 'sbc';
neuron(z).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-05-11-2/data014/elecResp_n797_p47_f5.mat';
neuron(z).scanElecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-05-11-2/data002/elecResp_n798_p47.mat';
neuron(z).PW = 100;
%for plotting ACF
neuron(z).WN.neuronsPath = '/snle/lab/Experiments/Array/Analysis/2011-05-11-2/data009/data009.neurons';
neuron(z).WN.id = 797;

z = z+1;
neuron(z).path = '/snle/lab/Experiments/Array/Analysis/2011-05-11-2/data015-filtered-cf-lh/data015-filtered-cf-lh.neurons';
neuron(z).movieTimePath = '/snle/lab/Experiments/Array/Analysis/2011-05-11-2/data015-filtered/data015-movieTimes.mat';
neuron(z).id = 110; %n141 in data009.ei
neuron(z).stimElec = 10;
neuron(z).type = 'offPar';
neuron(z).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-05-11-2/data015/elecResp_n141_p10_f5.mat';
neuron(z).scanElecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-05-11-2/data013/elecResp_n141_p10.mat';
neuron(z).PW = 100;
%for plotting ACF
neuron(z).WN.neuronsPath = '/snle/lab/Experiments/Array/Analysis/2011-05-11-2/data009/data009.neurons';
neuron(z).WN.id = 141;

% 30-micron array
% z = z+1;
% neuron(z).path = '/snle/lab/Experiments/Array/Analysis/2011-05-11-8/data007-filtered/data007-filtered.neurons';
% neuron(z).movieTimePath = '/snle/lab/Experiments/Array/Analysis/2011-05-11-8/data007-filtered/data007-movieTimes.mat';
% neuron(z).id = 573; %n556
% neuron(z).stimElec = 38;
% neuron(z).type = 'sbc';
% neuron(z).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-05-11-8/data007/elecResp_n556_p38_f5.mat';
% neuron(z).scanElecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-05-11-8/data005/elecResp_n556_p38.mat';
% neuron(z).PW = 100;

z = z+1;
neuron(z).path = '/snle/lab/Experiments/Array/Analysis/2011-06-24-0/data004-filtered/data004-filtered.neurons';
neuron(z).movieTimePath = '/snle/lab/Experiments/Array/Analysis/2011-06-24-0/data004-filtered/data004-movieTimes.mat';
neuron(z).id = 739; %n768 in data000.ei
neuron(z).stimElec = 52;
neuron(z).type = 'sbc';
neuron(z).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-06-24-0/data004/elecResp_n768_p52_f5.mat';
neuron(z).scanElecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-06-24-0/data002/elecResp_n768_p52.mat';
neuron(z).PW = 100;
%for plotting ACF
neuron(z).WN.neuronsPath = '/snle/lab/Experiments/Array/Analysis/2011-06-24-0/data000/data000.neurons';
neuron(z).WN.id = 768;

z = z+1;
neuron(z).path = '/snle/lab/Experiments/Array/Analysis/2011-06-24-5/data004-filtered-cf-lh/data004-filtered-cf-lh.neurons';
neuron(z).movieTimePath = '/snle/lab/Experiments/Array/Analysis/2011-06-24-5/data004-filtered/data004-movieTimes.mat';
neuron(z).id = 189; %n186 in data000.ei
neuron(z).stimElec = 10;
neuron(z).type = 'offMidg';
neuron(z).elecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-06-24-5/data004/elecResp_n186_p10_f5.mat';
neuron(z).scanElecRespPath = '/snle/lab/Experiments/Array/Analysis/2011-06-24-5/data002/elecResp_n186_p10.mat';
neuron(z).PW = 100;
%for plotting ACF
neuron(z).WN.neuronsPath = '/snle/lab/Experiments/Array/Analysis/2011-06-24-5/data000/data000.neurons';
neuron(z).WN.id = 186;


