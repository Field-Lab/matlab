function [xscale p] = xscalecr(cr, stim, varargin)
% XSCALECR  Calculate scaling necessary to match up C/R curves (Least Squares)
% usage: varargout = xscalecr(cr, stim, varargin)
%
% Unclear how this relates to NORMCDFXSCALE.  Maybe I wrote two different
% versions without realizing?  Or NORMCDFXSCALE was meant to replace? 
% Seems like NORMCDFXSCALE is the more developed at this point in any case.
%
% 2011 phli
%


opts = inputParser;
opts.addParamValue('fitfunc', @crnormcdffit);%@(p,x)(p(2)*normcdf(p(3).*abs(abs(x) - x) + p(4).*abs(x), 1, p(1))));
opts.addParamValue('styles', {'.' 'r.' 'k.'});
opts.addParamValue('MaxFunEvals', 5000);
opts.addParamValue('MaxIter', 5000);
opts.addParamValue('indices', 1:size(cr,3));
opts.addParamValue('plot', false);
opts.parse(varargin{:});
opts = opts.Results;


optimopts = optimset('MaxFunEvals', opts.MaxFunEvals, 'MaxIter', opts.MaxIter);
x = get_intensities(stim); 
fitx = (1.1*min(x)):.01:(1.1*max(x));


for ind = opts.indices;
    % First fit to summed trial
    trial = 3;
    y{trial} = cr(:,trial,ind);
    if max(y{trial}) == 0, continue; end
    p{ind} = lsqcurvefit(opts.fitfunc, [0.5 max(y{trial}) 1 5], x, y{trial}, [], [], optimopts);
    xscale{ind}(trial) = 1;
    
    trial = 1;
    y{trial} = cr(:,trial,ind);
    xscale{ind}(trial) = lsqcurvefit(@(xscale,x)(opts.fitfunc(p{ind},x.*xscale)), 1, x, y{trial});
    
    trial = 2;
    y{trial} = cr(:,trial,ind);
    xscale{ind}(trial) = lsqcurvefit(@(xscale,x)(opts.fitfunc(p{ind},x.*xscale)), 1, x, y{trial});
    

    % Plot that
    if opts.plot
        figure; hold on;

        fity = opts.fitfunc(p{ind},fitx);
        plot(fitx,fity,'k--');
        title(num2str(p{ind}));

        for trial = 1:3
            plot(x.*xscale{ind}(trial),y{trial},opts.styles{trial});
        end
    end
end



if nargout == 0 && opts.plot, clear xscale; end



function y = crnormcdffit(p,x)
xpos = x >= 0;
xneg = ~xpos;
y(xpos,1) = p(2) .* normcdf(p(3).*x(xpos),      1, p(1));
y(xneg,1) = p(2) .* normcdf(p(4).*abs(x(xneg)), 1, p(1));