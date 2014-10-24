function varargout = getOutputs(fn, ixsOutputs)

output_cell = cell(1,max(ixsOutputs));
[output_cell{:}] = (fn());
varargout = output_cell(ixsOutputs);