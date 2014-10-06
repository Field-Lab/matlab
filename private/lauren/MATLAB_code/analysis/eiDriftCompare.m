%% script to plot comparisons of EIs from different WN noise runs throughout an experiment
clear all
% 

% gain different between runs (? to 160)
dataPaths = {'2008-11-10-3/data000-PCA/data000-PCA','2008-11-10-3/data012/data012'};
savePath = '/snle/home/lhruby/Desktop/ei drift analysis/2008-11-10-3/';

pairsToCompare = [873 872;
                  274 287;
                  318 316;
                  648 722;
                  767 766;
                  873 872;
                  248 197;
                  392 333;
                  526 377;
                  768 767];

% dataPaths = {'2008-08-27-4/data000-NW-lh/data000-NW-lh','2008-08-27-4/data004-NW/data004-NW'};
% 
% pairsToCompare = [33 31;
%                   76 76;
%                   107 106;
%                   168 241;
%                   226 271;
%                   603 602;
%                   781 781;
%                   903 901;
%                   317 316;
%                   333 331;
%                   871 813;
%                   542 543];

% dataPaths = {'2008-08-26-0/data000-NW/data000-NW','2008-08-26-0/data007-NW/data007-NW'};
% 
% pairsToCompare = [109 137;
%                   227 227;
%                   557 559;
%                   889 828;
%                   198 197;
%                   725 635;
%                   81  1;
%                   261 259;
%                   333 333;
%                   346 317;
%                   406 406;
%                   631 631;
%                   347 347;
%                   301 301
%                   287 287];


% dataPaths = {'2009-09-05-1/data001/data001','2009-09-05-1/data009/data009'};
% 
% pairsToCompare = [948 901;
%                   196 166; 
%                   541 586;
%                   676 706;
%                   271 227;
%                   439 406;
%                   559 590];

% dataPaths = {'2009-09-03-1/data002/data002', '2009-09-03-1/data004/data004'};
% 
% pairsToCompare = [77  76;
%                   136 136;
%                   227 226;
%                   211 243;
%                   331 331;
%                   901 901
%                   1   2;
%                   138 92;
%                   829 831;
%                   286 334;
%                   168 168;
%                   287 333;
%                   347 393;
%                   453 499;
%                   586 588;
%                   722 676;
%                   902 904;
%                   272 272;
%                   663 663;
%                   766 603];

%dataPaths = {'2009-09-06-0/data001/data001', '2009-09-06-0/data010/data010'};
%
% pairsToCompare = [936 937;
%                   661 661;
%                   64  65;
%                   166 167;
%                   602 602;
%                   632 632;
%                   31  31;
%                   136 106;
%                   892 886;
%                   601 857;
%                   706 752;
%                   634 571];


% dataPaths = {'2010-03-05-3/data000/data000', '2010-03-05-3/data003/data003'};
% 
% pairsToCompare = [751 752;
%                  31 31
%                  468 423
%                  646 676
%                  661 661
%                  908 901
%                  182 227
%                  274 272];


% dataPaths = {'2010-03-05-4/data000-lh/data000-lh', '2010-03-05-4/data003-lh/data003-lh'};
% 
%off cells
% pairsToCompare = [872 871;
%                   602 213;
%                   861 856;
%                   828 20;
%                   786 811];
              
%on cells
% pairsToCompare = [917 931
%                   857 888
%                   722 721
%                   669 680
%                   557 481
%                   170 166];


nPairs = size(pairsToCompare, 1);
marks_params.strength = {'inner', [0 1 0]};

datarun1 = load_data(dataPaths{1});
datarun1 = load_sta(datarun1);
datarun1 = load_neurons(datarun1);
datarun1 = load_params(datarun1, 'verbose', true);
datarun1 = get_sta_summaries(datarun1, 'all', 'keep_rf_coms', false, 'marks_params', marks_params);


datarun2 = load_data(dataPaths{2});
datarun2 = load_sta(datarun2);
datarun2 = load_neurons(datarun2);
datarun2 = load_params(datarun2, 'verbose', true);
datarun2 = get_sta_summaries(datarun2, 'all', 'keep_rf_coms', false, 'marks_params', marks_params);


%%

for ii = 1:nPairs

    cell_id1 = pairsToCompare(ii, 1);
    cell_id2 = pairsToCompare(ii, 2);
    
    % receptive field comparison
    rfData1 = get_rf(datarun1,cell_id1,'polarity',true);
    rfData2 = get_rf(datarun2,cell_id2,'polarity',true);
    
    rfData1 = norm_image(rfData1);
    rfData2 = norm_image(rfData2);
    
    rfDiff = rfData2 - rfData1 + 0.5*ones(size(rfData1));
    rfDiff(rfDiff>1) = 1;
    rfDiff(rfDiff<0) = 0;
    
    figure('position', [100 100 600 350])
    
    axes('position', [0.05 0.45 0.25 0.38])
    image(rfData1)
    title(['WN run 1 (cell ' num2str(cell_id1) ')'])
    axis off
    
    axes('position', [0.35 0.45 0.25 0.38])
    image(rfData2)
    title(['WN run 2 (cell ' num2str(cell_id2) ')'])
    axis off
    text(0, 1.3, [dataPaths{1} 10 dataPaths{2}], 'units', 'normalized', 'FontSize', 8)
    
    axes('position', [0.65 0.45 0.25 0.38])
    image(rfDiff)
    title('difference')
    axis off
    
    % ei comparison    
    aEi1 = axes('position', [0.05 0.05 0.25 0.38]);
    [eiAmps1 scaleFactor] = plotEi61(datarun1.names.rrs_ei_path, cell_id1, 'axesH', aEi1);
    scaleFactor = scaleFactor*0.5;
    cla(aEi1)
    plotEi61(datarun1.names.rrs_ei_path, cell_id1, 'axesH', aEi1, 'manualScale', scaleFactor)
    
    aEi2 = axes('position', [0.35 0.05 0.25 0.38]);
    eiAmps2 = plotEi61(datarun2.names.rrs_ei_path, cell_id2, 'axesH', aEi2, 'manualScale', scaleFactor);
    
    eiAmpsDiffTemp = eiAmps2{1} - eiAmps1{1};
    eiAmpsDiff{1} = eiAmpsDiffTemp;
    eiAmpsDiff{1}(eiAmpsDiff{1}<0) = 0;
    eiAmpsDiff{2} = -eiAmpsDiffTemp;
    eiAmpsDiff{2}(eiAmpsDiff{2}<0) = 0;
    
    aEiDiff = axes('position', [0.65 0.05 0.25 0.38]);
    plotEi61('', 0, 'axesH', aEiDiff, 'manualScale', scaleFactor, 'eiData', eiAmpsDiff)
    
%     if exist('savePath', 'var')
%         saveas(gcf, [savePath 'cell' num2str(cell_id1) '.eps'], 'psc2')
%     end
end
