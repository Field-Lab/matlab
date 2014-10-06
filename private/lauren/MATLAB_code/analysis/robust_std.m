function rstd = robust_std(list)
% ROBUST_STD   compute robust standard deviation
%
% usage: rstd = robust_std(list)
%
% arguments:     list - vector of numbers
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


if 0

    %  LINE BY LINE TRANSLATION FROM LISP

    % (defmethod robust-standard-deviation ((list list))
    %   (let* ((copy (copy list))
    %          (median (median copy))
    %          median-absolute-deviation)
    %     (sub copy median :-> copy)
    list = list - median(list);

    %     (abs-value copy :-> copy)
    list = abs(list);

    %     (setq median-absolute-deviation (median copy))
    med_abs_dev = median(list);

    %     (/ median-absolute-deviation 0.6741891400433162)))
    rstd = med_abs_dev / 0.6741891400433162;

else
    %  MORE COMPACT VERSION
    rstd = median(abs(list - median(list))) / 0.6741891400433162;

end

