function datarun = calc_rf_snrs(datarun, varargin)
% CALC_RF_SNRS      SNR calculations on STA summary RFs
% usage: datarun = calc_rf_snrs(datarun, opts)
%
% Datarun wrapper for CALC_RF_SNRS; see low level for more documentation.
%
% See also: CALC_RF_SNRS
%
% 2012-09 phli
%

snrcalcs = calc_rf_snrs(datarun.stas.rfs, varargin{:});
datarun.stas.maxsigs = snrcalcs.maxsigs;
datarun.stas.numsigs = snrcalcs.numsigs;
datarun.stas.medsigs = snrcalcs.medsigs;