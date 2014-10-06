function rstd = robust_std(data)
% ROBUST_STD   compute robust standard deviation
%
% usage: rstd = robust_std(data)
%
% arguments:     data - vector of numbers 
%                       or array - operates along 1st dim 
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
%

if 0

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

else
    
    %  MORE COMPACT VERSION
    
    % if one dimensional, make sure it is a column
    if find(size(data)==1)==1
        data=data';
    end
    
    % identify median of each column, and subtract it from each entry in the column
    data = data - repmat(median(data),size(data,1),1);
    
    % rest of computation
    rstd = median(abs(data)) / 0.6741891400433162;

end

