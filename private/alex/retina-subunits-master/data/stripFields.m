function datarun = stripFields(datarun)

% remove the stas, and make spike counts uint8

datarun.stas = rmfield(datarun.stas,'stas');
datarun.spike_rate = uint8(datarun.spike_rate);
datarun.cone_inputs = single(datarun.cone_inputs);