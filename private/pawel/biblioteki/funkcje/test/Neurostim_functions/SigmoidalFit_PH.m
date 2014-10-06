function cf_=SigmoidalFit_PR(AllMovies,Effic)
%CREATEFIT    Create plot of datasets and fits
%   CREATEFIT(ALLMOVIES,EFFIC)
%   Creates a plot, similar to the plot in the main curve fitting
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with cftool
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  1
%   Number of fits:  1


% Data from dataset "Effic vs. AllMovies":
%    X = AllMovies:
%    Y = Effic:
%    Unweighted
%
% This function was automatically generated on 11-May-2010 16:38:11

% Set up figure to receive datasets and fits
figure(111);
f_ = clf;
figure(f_);
set(f_,'Units','Pixels','Position',[664 96 680 481]);
legh_ = []; legt_ = {};   % handles and text for legend
xlim_ = [Inf -Inf];       % limits of x axis
ax_ = axes;
set(ax_,'Units','normalized','OuterPosition',[0 0 1 1]);
set(ax_,'Box','on');
axes(ax_); hold on;


% --- Plot data originally in dataset "Effic vs. AllMovies"
AllMovies = AllMovies(:);
Effic = Effic(:);
h_ = line(AllMovies,Effic,'Parent',ax_,'Color',[0.333333 0 0.666667],...
    'LineStyle','none', 'LineWidth',1,...
    'Marker','.', 'MarkerSize',12);
xlim_(1) = min(xlim_(1),min(AllMovies));
xlim_(2) = max(xlim_(2),max(AllMovies));
legh_(end+1) = h_;
legt_{end+1} = 'Effic vs. AllMovies';

% Nudge axis limits beyond data limits
if all(isfinite(xlim_))
    xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
    set(ax_,'XLim',xlim_)
else
    set(ax_, 'XLim',[-0.5, 152.5]);
end


% --- Create fit "fit 3"
fo_ = fitoptions('method','NonlinearLeastSquares','Normalize','on');
ok_ = isfinite(AllMovies) & isfinite(Effic);
if ~all( ok_ )
    warning( 'GenerateMFile:IgnoringNansAndInfs', ...
        'Ignoring NaNs and Infs in data' );
end
st_ = [0.28946898665270315 0.5608943345498274 ];
set(fo_,'Startpoint',st_);
ft_ = fittype('100/(1+exp(-x*a+b))',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a', 'b'});

% Fit this model using new data
cf_ = fit(AllMovies(ok_),Effic(ok_),ft_,fo_)

% Or use coefficients from the original fit:
if 0
    cv_ = { 298.37809408436232, 357.97621432406191};
    cf_ = cfit(ft_,cv_{:});
end

% Plot this fit
h_ = plot(cf_,'fit',0.95);
legend off;  % turn off legend from plot method call
set(h_(1),'Color',[1 0 0],...
    'LineStyle','-', 'LineWidth',2,...
    'Marker','none', 'MarkerSize',6);
legh_(end+1) = h_(1);
legt_{end+1} = 'fit 3';

% Done plotting data and fits.  Now finish up loose ends.
hold off;
leginfo_ = {'Orientation', 'vertical', 'Location', 'NorthEast'};
h_ = legend(ax_,legh_,legt_,leginfo_{:});  % create legend
set(h_,'Interpreter','none');
xlabel(ax_,'');               % remove x label
ylabel(ax_,'');               % remove y label
