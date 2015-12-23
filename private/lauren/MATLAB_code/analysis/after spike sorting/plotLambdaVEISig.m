function plotLambdaVEISig(CI, varargin)


p = inputParser;

p.addRequired('CI', @isstruct)

p.addParamValue('excludeBadElecs', true, @islogical)
p.addParamValue('normEachCell', true, @islogical) %whether to normalize lambdas between 0 and 1 for each cell individually or for all cells together
p.addParamValue('allowNeg', true, @islogical) %whether to allow negative lambda values (vs. normalizing so that min(lambda) == 0)
p.addParamValue('useAbsEISig', true, @islogical)
p.addParamValue('excludeCellsInd', [], @isnumeric) %index (within cellInfo) of particular cell(s) to exclude from the analysis
p.addParamValue('nShuffleReps', 0, @isnumeric)
p.addParamValue('lambdaType', 'all', @ischar) %valid strings are 'all', 'positive', or 'negative': lambdas of unincluded sign are clipped to 0

p.parse(CI, varargin{:})

params = p.Results;


%% collect EI and lambda values for each cell
lam_all = {};
sElecEISig_all = {};

eiSig_data = [];
lam_data = [];

for ii = 1:length(CI)
    if ~any(params.excludeCellsInd == ii)
        sElecEISig = CI(ii).ei(CI(ii).details.sElecs,:);
        pElecEISig = CI(ii).ei(CI(ii).pElec,:);
        if params.useAbsEISig
            sElecEISig = max(abs(sElecEISig), [], 2);
            pElecEISig = max(abs(pElecEISig));
        else %use peak negative EI value
            sElecEISig = max(-sElecEISig, [], 2);
            pElecEISig = max(-pElecEISig);
        end
        
        %normalize to primary electrode signal
        sElecEISig = sElecEISig/pElecEISig;
        lam = CI(ii).details.lambdas;
        
        if params.excludeBadElecs
            goodElecs = ~(ismember(CI(ii).details.sElecs, CI(ii).badSElecs));
        else
            goodElecs = true(size(sElecEISig));
        end
        
        sElecEISig_all{end+1} = sElecEISig(goodElecs)'; %#ok<AGROW>
        lam_all{end+1} = lam(goodElecs); %#ok<AGROW>
        
        switch params.lambdaType
            case 'all' %leave lambda values as is
            case 'positive' %clip all negative values to 0
                lam_all{end}(lam_all{end}<0) = 0; %#ok<AGROW>
            case 'negative' %clip all positive values to 0
                lam_all{end}(lam_all{end}>0) = 0; %#ok<AGROW>
            otherwise
                error('invalid parameter value for ''lambdaType''')
        end
        
        %store actual data
        eiSig_data = [eiSig_data sElecEISig_all{end}];
        lam_data = [lam_data lam_all{end}];

       
        %plotting
        %figure; hold on
%         [~, iSort] = sort(sElecEISig);
%         plot(sElecEISig(iSort), lam(iSort), '-', 'color', lineColors(ii,:))
%         plot(sElecEISig(iSort), abs(lam(iSort)), '--', 'color', lineColors(ii,:))
%         for jj = 1:6
%             if ~any(CI(ii).badSElecs == CI(ii).details.sElecs(jj))
%                 plot(sElecEISig(jj), lam(jj), 'o', 'MarkerEdgeColor', lineColors(ii,:))
%                 plot(sElecEISig(jj), abs(lam(jj)), 'o', 'MarkerEdgeColor', lineColors(ii,:))
%             else
%                 plot(sElecEISig(jj), lam(jj), 'x', 'MarkerEdgeColor', [0 0 0])
%                 plot(sElecEISig(jj), abs(lam(jj)), 'x', 'MarkerEdgeColor', [0 0 0])
%             end
%         end
%         xPlotLim = get(gca, 'xlim');
%         yPlotLim = get(gca, 'ylim');
        
        %text(xPlotLim(1)+(xPlotLim(2)-xPlotLim(1))*0.1, yPlotLim(2)-(yPlotLim(2)-yPlotLim(1))*0.1)
        
%         plot(xPlotLim, [0 0], 'k--')
        
        %     title(['neuron ' num2str(CI(ii).id)])
        %     xlabel('EI amplitude, normalized to primary electrode')
        %     ylabel('lambda')
        %     keyboard
    end
end

%% shuffle analysis

r_sq_shuff = inf*ones(1,params.nShuffleReps);

xExt = [min(eiSig_data)  max(eiSig_data)];

figure
hold on
for ii = 1:params.nShuffleReps
    lam_shuff = [];
    eiSig_shuff = [];
    for jj = 1:length(lam_all)
        lam_shuff = [lam_shuff lam_all{jj}(randperm(length(lam_all{jj})))]; %#ok<AGROW> %randomly shuffle lambda values
        eiSig_shuff = [eiSig_shuff sElecEISig_all{jj}]; %#ok<AGROW>
    end
    
    %metric
    
    %linear regression
    A = ones(length(eiSig_shuff), 2);
    A(:,1) = eiSig_shuff;
    b = lam_shuff';
    x = A\b;
    m = x(1);
    y_0 = x(2);
    
    % calculate R-squared (where model is observed = predicted from additivity, vs. observed = predicted by best linear fit)
    pred = y_0 + m*eiSig_shuff;
    obs = lam_shuff;
    r_sq_shuff(ii) = 1 - sum((obs-pred).^2)/sum((obs-mean(obs)).^2);
    
    plot(eiSig_shuff, lam_shuff, '.', 'markerEdgeColor', [0.8 0.8 0.8])
    plot(xExt, xExt*m+y_0, '-', 'color', [0.8 0.8 0.8])
end



%compare to R^2 of data

A = ones(length(eiSig_data), 2);
A(:,1) = eiSig_data;
b = lam_data';
x = A\b;
m = x(1);
y_0 = x(2);

pred = y_0 + m*eiSig_data;
obs = lam_data;
r_sq_data = 1 - sum((obs-pred).^2)/sum((obs-mean(obs)).^2);

plot(eiSig_data, lam_data, 'r.')
plot(xExt, xExt*m+y_0, 'r-')

title(['mean R^2 value = ' num2str(mean(r_sq_shuff)), ', sig = ' num2str(std(r_sq_shuff)), 10,...
    'data R^2 value = ' num2str(r_sq_data)])


%% plot results of shuffle analysis to see if mean and sigma converged

%calculate cumulative means and standard deviations to look for convergence
cumMeanR_sq = zeros(1,params.nShuffleReps);
cumSigR_sq = zeros(1,params.nShuffleReps);

for ii = 1:params.nShuffleReps
    cumMeanR_sq(ii,:) = mean(r_sq_shuff(1:ii));
    cumSigR_sq(ii,:) = std(r_sq_shuff(1:ii));
end
figure
hold on
plot(cumMeanR_sq(:,ii), 'k-')
plot(cumSigR_sq(:,ii), 'k--')
xlabel('iteration')
ylabel('cumulative mean/standard deviation')
























