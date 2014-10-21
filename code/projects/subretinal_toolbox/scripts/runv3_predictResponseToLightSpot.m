clear;
savePlot = true;    % If true will save all plots created by predictResponseToLigt in figuresFolder

%% Specifying neuron IDs

largeRFsIDs = [212 273 301 452 467 526 542 571 646 796 857 871 991 ...
             1111 1126 1261 1262 1352 1398 1412 1426 1576 1592 1681 ...
             1831 1832 2027 2222 2374 2507 2716 2747 3079 3108 3451 ...
             3542 3886 3931 3993 4231 4322 4411 4486 4501 4636 4651 ...
             4726 4771 4921 5011 5056 5087 5147 5237 5297 5432 5776 ...
             5791 5806 5897 6001 6196 6211 6271 6436 6556 6616 6617 ...
             7007 7081 7201 7293 7366 7427 7636];
smallRFsIDs = [16 61 62 168 181 198 213 226 257 317 361 377 391 392 ...
             436 469 481 601 616 676 707 766 842 843 902 992 1007 ...
             1366 1501 1537 1561 1756 1816 1817 1891 1893 2088 2132 ...
             2252 2267 2312 2341 2656 2701 2778 2956 3185 3647 3811 ...
             3901 4086 4186 4203 5075 5192 5206 5536 5821 5851 5866 ...
             5956 5972 5986 6121 6227 6346 6466 6512 6526 6646 6691 ...
             6811 6813 6901 6902 6918 6931 6946 6991 7036 7082 7171 ...
             7202 7231 7291 7352 7411 7606];
allNeuronsIDs = [largeRFsIDs smallRFsIDs];

%% Moving the spot around a bit

dataFolder = '/media/MEA_PROCESSED_4/2012-08-06-0/data/data000/vision_processing/data000';
figuresFolder = '/media/MEA_PROCESSED_4/2012-08-06-0/data/data000/figures/predicted_response_allRFs';
spotCenterPosition = [-225 -70;
                      30 -70;
                      285 -70;
                      -225 110;
                      30 110;
                      -225 290;
                      30 290;
                      285 290];
                  
spotCenterPosition = [-45*ones(1,11) 255*ones(1,11);
                      linspace(-330,390,11) linspace(-330,390,11)].';
                  
allSpotDiameters = [10 20 40 70 100 140 210 280 350 420 490 560 700 840];
% allSpotDiameters = 140;

npoints = 5000;
maxdist = 2000;
f = @(p,x) p(1).*exp(-abs(x.^p(3))./p(2)^2);
f_inv = @(p,y) (-p(2)^2.*log(y./p(1))).^(1/p(3));
x = linspace(0,maxdist,npoints);

allThresholds = zeros(size(spotCenterPosition,1),length(allSpotDiameters));
for kk=1:size(spotCenterPosition)
    for ll=1:length(allSpotDiameters)
        spotDiameter = allSpotDiameters(ll);
        p = predictResponseToLightSpot(dataFolder, figuresFolder, spotCenterPosition(kk,:), ...
            spotDiameter,'neuronList',allNeuronsIDs,'savePlot',savePlot);
        y = f(p,x);
        allThresholds(kk,ll) = x(find(y<y(1)*0.5,1,'first'));
    end
end

%% With a select spot, change RF diameter, get 50% threshold and plot how 
% threshold changes vs spot size

meanThreshold = mean(allThresholds,1);

fh = figure(2); clf; set(fh,'color','white');
hold on
% plot([0 allSpotDiameters], [0 allSpotDiameters]*line(1) + line(2), '--b','linewidth',2)
plot(allSpotDiameters,meanThreshold,'r','linewidth',2);
for kk=1:size(allThresholds,1)
    scatter(allSpotDiameters,allThresholds(kk,:),'.k');
end

xlabel('Spot diameter, \mum')
ylabel('50% stimulation spread, \mum')
% legend('Linear fit of the spread','Mean stimulation spread','Spread for one spot size position')
legend('Mean stimulation spread, all RF sizes','Spread for one spot size position')
title('Stimulation spread vs. spot diameter, estimated from LE rat RF size','fontsize',12)
axis([0 max(allSpotDiameters) 0 1.1*max(allThresholds(:))])