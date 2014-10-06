function varargout = rf_snr_classification(datarun)
% RF_SNR_CLASSIFICATION     Run SNR calcs on RFs and spit out classification.txt
% usage: [classifications] = rf_snr_classification(datarun)
%
% Currently not very flexible.  Runs on all RFs it finds.  Writes three
% sets of decile classifications into the rrs-prefix folder.
%
% See also: WRITE_TXT_CLASSIFICATION, RESHAPE_TXT_CLASSIFICATION,
% CALC_RF_SNRS
%
% 2011-05 phli
%

if ~isfield(datarun, 'stas') || ~isfield(datarun.stas, 'rfs')
    warning('No RFs in datarun.  Run: load_sta(datarun, ''save_rf'', true)');
    return;
end


snr_calcs = calc_rf_snrs(datarun.stas.rfs);

maxclass = write_decile_classification('maxsigs', snr_calcs.maxsigs, datarun);
medclass = write_decile_classification('medsigs', snr_calcs.medsigs, datarun);
numclass = write_decile_classification('numsigs', snr_calcs.numsigs, datarun);

if nargout > 0
    varargout{1} = struct();
    varargout{1}.maxsigsclass = maxclass;
    varargout{1}.medsigsclass = medclass;
    varargout{1}.numsigsclass = numclass;
end


% Sort values into decile classes and write output
function classification = write_decile_classification(name, vals, datarun)
vals(isnan(vals)) = 0; % Otherwise NaN sorts to the top

classification = decile_class(name);

[~, sortind] = sort(vals);
for i = 1:length(sortind)
    cell_id = datarun.cell_ids(sortind(i));
    dec = ceil(i / length(sortind) * 10);
    classification.subclasses(dec).cells = [classification.subclasses(dec).cells cell_id];
end

write_txt_classification(classification, [datarun.names.rrs_prefix '_' name '_classification.txt']);



% Init an empty decile classification
function classification = decile_class(name)
classification = struct('name', 'All', 'cells', []);
for i = 1:10
    classification.subclasses(i).name  = [name num2str(i)];
    classification.subclasses(i).cells = [];
    classification.subclasses(i).subclasses = struct([]);
end