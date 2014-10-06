function snr_calcs = calc_rf_snrs(rfs, varargin)
% CALC_RF_SNRS      SNR calculations on STA summary RFs
% usage: snr_calcs = calc_rf_snrs(rfs, opts)
%
% RFs is a cell array of STA summary RFs; must currently be BW RFs.
%
% Outputs struct with three vector fields:
%   MAXSIGS, the maximum sigma value for each RF
%   MEDSIGS, the median sigma value within significant pixels for each RF
%   NUMSIGS, the number of significant pixels for each RF
%
% Significant pixels are those whose sigma value (calculated with
% ROBUST_STD) are above a cutoff threshold.
%
% OPTS:
%   robust_std_method   3   See ROBUST_STD for different methods
%
%   fieldsize           []  Used with FALSEPOS to determine the sigma
%                           threshold cutoff for significant pixels. If
%                           left blank, derived from numel(rfs{1}).
%
%   falsepos            1   Number of false positive significant stixels to
%                           accept in setting the sigma threshold. Used in
%                           conjunction with FIELDSIZE.
%
%   pthresh             []  If set, overrides FIELDSIZE and FALSEPOS to
%                           directly set the probability threshold used to
%                           calculate the sigma threshold.
%
%   disp                false
%   
%
% See also: LOAD_STA, RF_FROM_STA, ROBUST_STD
%
% 2011-05 phli
%

opts = inputParser;
opts.addParamValue('robust_std_method', 5);
opts.addParamValue('fieldsize', []);
opts.addParamValue('falsepos', 1);
opts.addParamValue('pthresh', []);
opts.addParamValue('sigthresh', []);
opts.addParamValue('disp', false);
opts.parse(varargin{:});
opts = opts.Results;


% Determine the sigma threshold to use
if isempty(opts.sigthresh)
    if isempty(opts.pthresh)
        if isempty(opts.fieldsize)
            opts.fieldsize = numel(rfs{1});
        end
        opts.pthresh = 1 - (opts.falsepos / opts.fieldsize);
    end
    opts.sigthresh = norminv(opts.pthresh, 0, 1);
end


% Do RF calculations
numrfs = length(rfs);
numsigs = nan(1,numrfs);
maxsigs = nan(1,numrfs);
medsigs = nan(1,numrfs);
for i = 1:numrfs
    rf = rfs{i};
    if isempty(rf), continue; end
    
    rstd = robust_std(rf(:), opts.robust_std_method);
    normrf = rf ./ rstd;
    maxsigs(i) = max(normrf(:));

    sigpix = normrf > opts.sigthresh;
    sigvals = normrf(sigpix);
    numsigs(i) = length(sigvals);    
    if numsigs(i) > 0, medsigs(i) = fast_median(sigvals); end
    
    if opts.disp, imagesc(sigpix); colormap gray; end
end


% Compose output struct
snr_calcs.maxsigs = maxsigs;
snr_calcs.medsigs = medsigs;
snr_calcs.numsigs = numsigs;