function rstd = robust_std(data, method)
% ROBUST_STD   compute robust standard deviation
%
% usage: rstd = robust_std(data, [method])
%
% arguments:     data - vector of numbers 
%                       or array - operates along 1st dim
%
%                method - method for doing calculation. In general only
%                         method 1 (old method using Matlab median) or
%                         method 3 (new faster method requiring compiled
%                         fast_median.mex file) are used. See code for more
%                         details.
%
% outputs:       rstd - robust standard deviation
%
%
% translated from lisp (statistics.lisp) to matlab by JLG, 2007-07-17
%
% ;; Function ROBUST-STANDARD-DEVIATION
% ;; ESF 2003-08-13
% ;; made more efficient by EJC 2003-08-14
% ;; calculates a robust standard deviation
% ;; see PJ Huber (1981), Robust Statistics, Wiley, call number QA 276.H785
% ;; hardcoded constant is inverse cumulative normal evaluated at 0.75
% ;; this computation performs well if
% ;; (1) data near the mean approximately follow a normal distribution,
% ;; (2) the outliers are outside the 25%-75% percentiles.
%
% 2008-10 greschner loop over 2nd dim
% 2008-10 gauthier  remove loop for increased speed
% 2010-11 phli      added additional methods, ability to choose method, default to old method
%

if nargin < 2
    method = 1;
end

switch method
    case -1
        %%% DOESN'T WORK; SEEMS TO BE MISSING ITS LOOP!
        
        %  LINE BY LINE TRANSLATION FROM LISP  
        
        % only works when data is a vector
        
        % (defmethod robust-standard-deviation ((list list))
        %   (let* ((copy (copy list))
        %          (median (median copy))
        %          median-absolute-deviation)
        %     (sub copy median :-> copy)
        data = data - median(data);
        
        %     (abs-value copy :-> copy)
        data = abs(data);
        
        %     (setq median-absolute-deviation (median copy))
        med_abs_dev = median(data);
        
        %     (/ median-absolute-deviation 0.6741891400433162)))
        rstd = med_abs_dev / 0.6741891400433162;
        
        
        
    case 1
        
        %  MORE COMPACT VERSION
        
        % if one dimensional, make sure it is a column
        if find(size(data)==1)==1
            data=data';
        end
        
        % identify median of each column, and subtract it from each entry in the column
        data = data - repmat(median(data),size(data,1),1);
        
        % rest of computation
        rstd = median(abs(data)) ./ 0.6741891400433162;
        
    case 2
        % Use the statistics toolbox MAD() function.
        % Unfortunately, this is not particularly fast either
        
        % if one dimensional, make sure it is a column
        if find(size(data)==1)==1
            data=data';
        end
        
        rstd = mad(data, 1) ./ 0.6741891400433162;
        
    case 3
        % Use fast_median, derived from C++ nth_element, which uses
        % quickselect to get median, faster than MathWorks median which
        % simply sorts the whole thing!
        
        % if one dimensional, make sure it is a column
        if find(size(data)==1)==1
            data=data';
        end
        
        % identify median of each column, and subtract it from each entry in the column
        data = data - repmat(fast_median(data),size(data,1),1);
        
        % rest of computation
        rstd = fast_median(abs(data)) ./ 0.6741891400433162;
    
    case 4
        % MatLab native optimized version; avoids calling median twice and therefore
        % sorting twice.  Instead calls sort once, gets median as median
        % function normally does, then gets median of absolute values using
        % the already semi-sorted data.
        %
        % This is faster for lists up to around 10,000 long, but for longer
        % it bogs down, probably due to MatLab loops; should try moving
        % sorted list merge operation into C.
        
        % if one dimensional, make sure it is a column
        if find(size(data)==1)==1
            data=data';
        end
        
        % Sort ONCE
        data = sort(data, 1);
        
        % Get the first median, essentially same algorithm as used by
        % MathWorks MEDIAN()
        nCompare = size(data, 1);
        half = floor(nCompare/2);
        meds = data(half+1, :);
        if 2*half == nCompare        % Average if even number of elements
            meds = data(half, :) .* 0.5 + meds .* 0.5;
        end

        % Subtract the median, take absolute value
        data = abs(data - repmat(meds, size(data, 1), 1));

        
        % Instead of sorting again to get the median this time, use the
        % fact that the data on either side of the median value are already
        % sorted.
        
        %rstd = zeros(1, size(data, 2));

        % Treat the values before the median and after the median as two 
        % separate, already sorted, positive arrays.  Run through them
        % pulling off the mins in order until we get to the new median.
        for i = 1:size(data, 2)
            
            % Run through the mins until we get to two before the median
            comps = 0;
            posind = half+1;
            negind = half;
            while comps < half - 1
                if data(negind, i) < data(posind, i)
                    negind = negind - 1;
                else
                    posind = posind + 1;
                end
                
                comps = comps + 1;
            end

            % If even, get the next two mins and average them, else just
            % get the min two ahead.
            if 2*half == nCompare
                if data(negind, i) < data(posind, i)
                    first = data(negind, i);
                    negind = negind - 1;
                else
                    first = data(posind, i);
                    posind = posind + 1;
                end
                
                if data(negind, i) < data(posind, i)
                    second = data(negind, i);
                else
                    second = data(posind, i);
                end
                
                med = first .* 0.5 + second .* 0.5;
            else
                if data(negind, i) < data(posind, i)
                    negind = negind - 1;
                else
                    posind = posind + 1;
                end
                
                if data(negind, i) < data(posind, i)
                    med = data(negind, i);
                else
                    med = data(posind, i);
                end
            end
            
            rstd(i) = med;
        end
        
        rstd = rstd ./ 0.6741891400433162;
        
        
        
    case {5 6}
        % This version is even faster than 3, but it edits the input matrix
        % in place, which could be dangerous.  Initial testing suggests we
        % have accounted for the danger though.
        
        % Parallelism is now included in the normal functions if compiled with
        % OpenMP support, so no separate option 6.
        
        % if one dimensional, make sure it is a column
        if find(size(data)==1)==1
            data=data';
        end
        
        % identify median of each column, and subtract it from each entry in the column
        data = data - repmat(fast_median_ip(data),size(data,1),1);
        
        % rest of computation
        rstd = fast_median_ip(abs(data)) ./ 0.6741891400433162;
end
