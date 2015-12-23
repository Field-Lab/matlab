function datarun = load_electrode_position(datarun, varargin)
% load_electrode_position    fills position field, without load_ei
%
%              params - array_type = 61 | 512 | 519  default 512, uses fake id
%
% outputs:    datarun - datarun struct with the following fields added, as possible
%
%       datarun.positions
%
% greschner 


% SET UP OPTIONAL ARGUMENTS
    p = inputParser;
    p.addParamValue('array_type', 512);% 
    % parse inputs
    p.parse(varargin{:});
    params = p.Results;

if params.array_type==61
    temp=[9.0000000e+01   1.8000000e+02   1.2000000e+02   3.0000000e+01   1.5000000e+02   6.0000000e+01   9.0000000e+01   1.2000000e+02  -2.4000000e+02   6.0000000e+01   3.0000000e+01   0.0000000e+00   0.0000000e+00  -3.0000000e+01  -3.0000000e+01  -6.0000000e+01  -6.0000000e+01  -9.0000000e+01  -1.2000000e+02  -1.5000000e+02  -1.2000000e+02  -9.0000000e+01  -1.8000000e+02  -1.5000000e+02   2.4000000e+02  -6.0000000e+01  -2.1000000e+02  -1.2000000e+02  -1.8000000e+02  -2.4000000e+02  -2.1000000e+02  -1.5000000e+02  -9.0000000e+01  -1.8000000e+02  -1.2000000e+02  -3.0000000e+01  -1.5000000e+02  -6.0000000e+01  -9.0000000e+01  -1.2000000e+02   0.0000000e+00  -6.0000000e+01  -3.0000000e+01   0.0000000e+00   0.0000000e+00   3.0000000e+01   3.0000000e+01   6.0000000e+01   6.0000000e+01   9.0000000e+01   1.2000000e+02   1.5000000e+02   1.2000000e+02   9.0000000e+01   1.8000000e+02   1.5000000e+02   2.4000000e+02   6.0000000e+01   2.1000000e+02   1.2000000e+02   1.8000000e+02   2.4000000e+02   2.1000000e+02   1.5000000e+02;...
   6.0000000e+01   1.2000000e+02   1.2000000e+02   6.0000000e+01   1.8000000e+02   1.2000000e+02   1.8000000e+02   2.4000000e+02   2.4000000e+02   2.4000000e+02   1.8000000e+02   1.2000000e+02   2.4000000e+02   1.8000000e+02   6.0000000e+01   2.4000000e+02   1.2000000e+02   1.8000000e+02   2.4000000e+02   1.8000000e+02   1.2000000e+02   6.0000000e+01   1.2000000e+02   6.0000000e+01   2.4000000e+02   0.0000000e+00   6.0000000e+01   0.0000000e+00   0.0000000e+00   0.0000000e+00  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -1.2000000e+02  -1.2000000e+02  -6.0000000e+01  -1.8000000e+02  -1.2000000e+02  -1.8000000e+02  -2.4000000e+02   0.0000000e+00  -2.4000000e+02  -1.8000000e+02  -1.2000000e+02  -2.4000000e+02  -1.8000000e+02  -6.0000000e+01  -2.4000000e+02  -1.2000000e+02  -1.8000000e+02  -2.4000000e+02  -1.8000000e+02  -1.2000000e+02  -6.0000000e+01  -1.2000000e+02  -6.0000000e+01  -2.4000000e+02   0.0000000e+00  -6.0000000e+01   0.0000000e+00   0.0000000e+00   0.0000000e+00   6.0000000e+01   6.0000000e+01];
end
    
if params.array_type==512
    temp=[4.6500000e+02   4.6500000e+02   4.6500000e+02   4.6500000e+02   4.3500000e+02   4.3500000e+02   4.3500000e+02   4.3500000e+02   4.0500000e+02   4.0500000e+02   4.0500000e+02   4.0500000e+02   3.7500000e+02   3.7500000e+02   3.7500000e+02   3.7500000e+02   3.4500000e+02   3.4500000e+02   3.4500000e+02   3.4500000e+02   3.1500000e+02   3.1500000e+02   3.1500000e+02   3.1500000e+02   2.8500000e+02   2.8500000e+02   2.8500000e+02   2.8500000e+02   2.5500000e+02   2.5500000e+02   2.5500000e+02   2.5500000e+02   2.2500000e+02   2.2500000e+02   2.2500000e+02   2.2500000e+02   1.9500000e+02   1.9500000e+02   1.9500000e+02   1.9500000e+02   1.6500000e+02   1.6500000e+02   1.6500000e+02   1.6500000e+02   1.3500000e+02   1.3500000e+02   1.3500000e+02   1.3500000e+02   1.0500000e+02   1.0500000e+02   1.0500000e+02   1.0500000e+02   7.5000000e+01   7.5000000e+01   7.5000000e+01   7.5000000e+01   4.5000000e+01   4.5000000e+01   4.5000000e+01   4.5000000e+01   1.5000000e+01   1.5000000e+01   1.5000000e+01   1.5000000e+01  -1.5000000e+01  -1.5000000e+01  -1.5000000e+01  -1.5000000e+01  -4.5000000e+01  -4.5000000e+01  -4.5000000e+01  -4.5000000e+01  -7.5000000e+01  -7.5000000e+01  -7.5000000e+01  -7.5000000e+01  -1.0500000e+02  -1.0500000e+02  -1.0500000e+02  -1.0500000e+02  -1.3500000e+02  -1.3500000e+02  -1.3500000e+02  -1.3500000e+02  -1.6500000e+02  -1.6500000e+02  -1.6500000e+02  -1.6500000e+02  -1.9500000e+02  -1.9500000e+02  -1.9500000e+02  -1.9500000e+02  -2.2500000e+02  -2.2500000e+02  -2.2500000e+02  -2.2500000e+02  -2.5500000e+02  -2.5500000e+02  -2.5500000e+02  -2.5500000e+02  -2.8500000e+02  -2.8500000e+02  -2.8500000e+02  -2.8500000e+02  -3.1500000e+02  -3.1500000e+02  -3.1500000e+02  -3.1500000e+02  -3.4500000e+02  -3.4500000e+02  -3.4500000e+02  -3.4500000e+02  -3.7500000e+02  -3.7500000e+02  -3.7500000e+02  -3.7500000e+02  -4.0500000e+02  -4.0500000e+02  -4.0500000e+02  -4.0500000e+02  -4.3500000e+02  -4.3500000e+02  -4.3500000e+02  -4.3500000e+02  -4.6500000e+02  -4.6500000e+02  -4.6500000e+02  -4.6500000e+02  -9.4500000e+02  -8.8500000e+02  -8.2500000e+02  -7.6500000e+02  -7.0500000e+02  -6.4500000e+02  -5.8500000e+02  -5.2500000e+02  -9.1500000e+02  -8.5500000e+02  -7.9500000e+02  -7.3500000e+02  -6.7500000e+02  -6.1500000e+02  -5.5500000e+02  -4.9500000e+02  -9.4500000e+02  -8.8500000e+02  -8.2500000e+02  -7.6500000e+02  -7.0500000e+02  -6.4500000e+02  -5.8500000e+02  -5.2500000e+02  -9.1500000e+02  -8.5500000e+02  -7.9500000e+02  -7.3500000e+02  -6.7500000e+02  -6.1500000e+02  -5.5500000e+02  -4.9500000e+02  -9.4500000e+02  -8.8500000e+02  -8.2500000e+02  -7.6500000e+02  -7.0500000e+02  -6.4500000e+02  -5.8500000e+02  -5.2500000e+02  -9.1500000e+02  -8.5500000e+02  -7.9500000e+02  -7.3500000e+02  -6.7500000e+02  -6.1500000e+02  -5.5500000e+02  -4.9500000e+02  -9.4500000e+02  -8.8500000e+02  -8.2500000e+02  -7.6500000e+02  -7.0500000e+02  -6.4500000e+02  -5.8500000e+02  -5.2500000e+02  -9.1500000e+02  -8.5500000e+02  -7.9500000e+02  -7.3500000e+02  -6.7500000e+02  -6.1500000e+02  -5.5500000e+02  -4.9500000e+02  -9.4500000e+02  -8.8500000e+02  -8.2500000e+02  -7.6500000e+02  -7.0500000e+02  -6.4500000e+02  -5.8500000e+02  -5.2500000e+02  -9.1500000e+02  -8.5500000e+02  -7.9500000e+02  -7.3500000e+02  -6.7500000e+02  -6.1500000e+02  -5.5500000e+02  -4.9500000e+02  -9.4500000e+02  -8.8500000e+02  -8.2500000e+02  -7.6500000e+02  -7.0500000e+02  -6.4500000e+02  -5.8500000e+02  -5.2500000e+02  -9.1500000e+02  -8.5500000e+02  -7.9500000e+02  -7.3500000e+02  -6.7500000e+02  -6.1500000e+02  -5.5500000e+02  -4.9500000e+02  -9.4500000e+02  -8.8500000e+02  -8.2500000e+02  -7.6500000e+02  -7.0500000e+02  -6.4500000e+02  -5.8500000e+02  -5.2500000e+02  -9.1500000e+02  -8.5500000e+02  -7.9500000e+02  -7.3500000e+02  -6.7500000e+02  -6.1500000e+02  -5.5500000e+02  -4.9500000e+02  -9.4500000e+02  -8.8500000e+02  -8.2500000e+02  -7.6500000e+02  -7.0500000e+02  -6.4500000e+02  -5.8500000e+02  -5.2500000e+02  -9.1500000e+02  -8.5500000e+02  -7.9500000e+02  -7.3500000e+02  -6.7500000e+02  -6.1500000e+02  -5.5500000e+02  -4.9500000e+02  -4.6500000e+02  -4.6500000e+02  -4.6500000e+02  -4.6500000e+02  -4.3500000e+02  -4.3500000e+02  -4.3500000e+02  -4.3500000e+02  -4.0500000e+02  -4.0500000e+02  -4.0500000e+02  -4.0500000e+02  -3.7500000e+02  -3.7500000e+02  -3.7500000e+02  -3.7500000e+02  -3.4500000e+02  -3.4500000e+02  -3.4500000e+02  -3.4500000e+02  -3.1500000e+02  -3.1500000e+02  -3.1500000e+02  -3.1500000e+02  -2.8500000e+02  -2.8500000e+02  -2.8500000e+02  -2.8500000e+02  -2.5500000e+02  -2.5500000e+02  -2.5500000e+02  -2.5500000e+02  -2.2500000e+02  -2.2500000e+02  -2.2500000e+02  -2.2500000e+02  -1.9500000e+02  -1.9500000e+02  -1.9500000e+02  -1.9500000e+02  -1.6500000e+02  -1.6500000e+02  -1.6500000e+02  -1.6500000e+02  -1.3500000e+02  -1.3500000e+02  -1.3500000e+02  -1.3500000e+02  -1.0500000e+02  -1.0500000e+02  -1.0500000e+02  -1.0500000e+02  -7.5000000e+01  -7.5000000e+01  -7.5000000e+01  -7.5000000e+01  -4.5000000e+01  -4.5000000e+01  -4.5000000e+01  -4.5000000e+01  -1.5000000e+01  -1.5000000e+01  -1.5000000e+01  -1.5000000e+01   1.5000000e+01   1.5000000e+01   1.5000000e+01   1.5000000e+01   4.5000000e+01   4.5000000e+01   4.5000000e+01   4.5000000e+01   7.5000000e+01   7.5000000e+01   7.5000000e+01   7.5000000e+01   1.0500000e+02   1.0500000e+02   1.0500000e+02   1.0500000e+02   1.3500000e+02   1.3500000e+02   1.3500000e+02   1.3500000e+02   1.6500000e+02   1.6500000e+02   1.6500000e+02   1.6500000e+02   1.9500000e+02   1.9500000e+02   1.9500000e+02   1.9500000e+02   2.2500000e+02   2.2500000e+02   2.2500000e+02   2.2500000e+02   2.5500000e+02   2.5500000e+02   2.5500000e+02   2.5500000e+02   2.8500000e+02   2.8500000e+02   2.8500000e+02   2.8500000e+02   3.1500000e+02   3.1500000e+02   3.1500000e+02   3.1500000e+02   3.4500000e+02   3.4500000e+02   3.4500000e+02   3.4500000e+02   3.7500000e+02   3.7500000e+02   3.7500000e+02   3.7500000e+02   4.0500000e+02   4.0500000e+02   4.0500000e+02   4.0500000e+02   4.3500000e+02   4.3500000e+02   4.3500000e+02   4.3500000e+02   4.6500000e+02   4.6500000e+02   4.6500000e+02   4.6500000e+02   5.2500000e+02   5.8500000e+02   6.4500000e+02   7.0500000e+02   7.6500000e+02   8.2500000e+02   8.8500000e+02   9.4500000e+02   4.9500000e+02   5.5500000e+02   6.1500000e+02   6.7500000e+02   7.3500000e+02   7.9500000e+02   8.5500000e+02   9.1500000e+02   5.2500000e+02   5.8500000e+02   6.4500000e+02   7.0500000e+02   7.6500000e+02   8.2500000e+02   8.8500000e+02   9.4500000e+02   4.9500000e+02   5.5500000e+02   6.1500000e+02   6.7500000e+02   7.3500000e+02   7.9500000e+02   8.5500000e+02   9.1500000e+02   5.2500000e+02   5.8500000e+02   6.4500000e+02   7.0500000e+02   7.6500000e+02   8.2500000e+02   8.8500000e+02   9.4500000e+02   4.9500000e+02   5.5500000e+02   6.1500000e+02   6.7500000e+02   7.3500000e+02   7.9500000e+02   8.5500000e+02   9.1500000e+02   5.2500000e+02   5.8500000e+02   6.4500000e+02   7.0500000e+02   7.6500000e+02   8.2500000e+02   8.8500000e+02   9.4500000e+02   4.9500000e+02   5.5500000e+02   6.1500000e+02   6.7500000e+02   7.3500000e+02   7.9500000e+02   8.5500000e+02   9.1500000e+02   5.2500000e+02   5.8500000e+02   6.4500000e+02   7.0500000e+02   7.6500000e+02   8.2500000e+02   8.8500000e+02   9.4500000e+02   4.9500000e+02   5.5500000e+02   6.1500000e+02   6.7500000e+02   7.3500000e+02   7.9500000e+02   8.5500000e+02   9.1500000e+02   5.2500000e+02   5.8500000e+02   6.4500000e+02   7.0500000e+02   7.6500000e+02   8.2500000e+02   8.8500000e+02   9.4500000e+02   4.9500000e+02   5.5500000e+02   6.1500000e+02   6.7500000e+02   7.3500000e+02   7.9500000e+02   8.5500000e+02   9.1500000e+02   5.2500000e+02   5.8500000e+02   6.4500000e+02   7.0500000e+02   7.6500000e+02   8.2500000e+02   8.8500000e+02   9.4500000e+02   4.9500000e+02   5.5500000e+02   6.1500000e+02   6.7500000e+02   7.3500000e+02   7.9500000e+02   8.5500000e+02   9.1500000e+02   5.2500000e+02   5.8500000e+02   6.4500000e+02   7.0500000e+02   7.6500000e+02   8.2500000e+02   8.8500000e+02   9.4500000e+02   4.9500000e+02   5.5500000e+02   6.1500000e+02   6.7500000e+02   7.3500000e+02   7.9500000e+02   8.5500000e+02   9.1500000e+02;...
        -3.0000000e+01  -1.5000000e+02  -2.7000000e+02  -3.9000000e+02  -9.0000000e+01  -2.1000000e+02  -3.3000000e+02  -4.5000000e+02  -3.0000000e+01  -1.5000000e+02  -2.7000000e+02  -3.9000000e+02  -9.0000000e+01  -2.1000000e+02  -3.3000000e+02  -4.5000000e+02  -3.0000000e+01  -1.5000000e+02  -2.7000000e+02  -3.9000000e+02  -9.0000000e+01  -2.1000000e+02  -3.3000000e+02  -4.5000000e+02  -3.0000000e+01  -1.5000000e+02  -2.7000000e+02  -3.9000000e+02  -9.0000000e+01  -2.1000000e+02  -3.3000000e+02  -4.5000000e+02  -3.0000000e+01  -1.5000000e+02  -2.7000000e+02  -3.9000000e+02  -9.0000000e+01  -2.1000000e+02  -3.3000000e+02  -4.5000000e+02  -3.0000000e+01  -1.5000000e+02  -2.7000000e+02  -3.9000000e+02  -9.0000000e+01  -2.1000000e+02  -3.3000000e+02  -4.5000000e+02  -3.0000000e+01  -1.5000000e+02  -2.7000000e+02  -3.9000000e+02  -9.0000000e+01  -2.1000000e+02  -3.3000000e+02  -4.5000000e+02  -3.0000000e+01  -1.5000000e+02  -2.7000000e+02  -3.9000000e+02  -9.0000000e+01  -2.1000000e+02  -3.3000000e+02  -4.5000000e+02  -3.0000000e+01  -1.5000000e+02  -2.7000000e+02  -3.9000000e+02  -9.0000000e+01  -2.1000000e+02  -3.3000000e+02  -4.5000000e+02  -3.0000000e+01  -1.5000000e+02  -2.7000000e+02  -3.9000000e+02  -9.0000000e+01  -2.1000000e+02  -3.3000000e+02  -4.5000000e+02  -3.0000000e+01  -1.5000000e+02  -2.7000000e+02  -3.9000000e+02  -9.0000000e+01  -2.1000000e+02  -3.3000000e+02  -4.5000000e+02  -3.0000000e+01  -1.5000000e+02  -2.7000000e+02  -3.9000000e+02  -9.0000000e+01  -2.1000000e+02  -3.3000000e+02  -4.5000000e+02  -3.0000000e+01  -1.5000000e+02  -2.7000000e+02  -3.9000000e+02  -9.0000000e+01  -2.1000000e+02  -3.3000000e+02  -4.5000000e+02  -3.0000000e+01  -1.5000000e+02  -2.7000000e+02  -3.9000000e+02  -9.0000000e+01  -2.1000000e+02  -3.3000000e+02  -4.5000000e+02  -3.0000000e+01  -1.5000000e+02  -2.7000000e+02  -3.9000000e+02  -9.0000000e+01  -2.1000000e+02  -3.3000000e+02  -4.5000000e+02  -3.0000000e+01  -1.5000000e+02  -2.7000000e+02  -3.9000000e+02  -9.0000000e+01  -2.1000000e+02  -3.3000000e+02  -4.5000000e+02  -4.5000000e+02  -4.5000000e+02  -4.5000000e+02  -4.5000000e+02  -4.5000000e+02  -4.5000000e+02  -4.5000000e+02  -4.5000000e+02  -3.9000000e+02  -3.9000000e+02  -3.9000000e+02  -3.9000000e+02  -3.9000000e+02  -3.9000000e+02  -3.9000000e+02  -3.9000000e+02  -3.3000000e+02  -3.3000000e+02  -3.3000000e+02  -3.3000000e+02  -3.3000000e+02  -3.3000000e+02  -3.3000000e+02  -3.3000000e+02  -2.7000000e+02  -2.7000000e+02  -2.7000000e+02  -2.7000000e+02  -2.7000000e+02  -2.7000000e+02  -2.7000000e+02  -2.7000000e+02  -2.1000000e+02  -2.1000000e+02  -2.1000000e+02  -2.1000000e+02  -2.1000000e+02  -2.1000000e+02  -2.1000000e+02  -2.1000000e+02  -1.5000000e+02  -1.5000000e+02  -1.5000000e+02  -1.5000000e+02  -1.5000000e+02  -1.5000000e+02  -1.5000000e+02  -1.5000000e+02  -9.0000000e+01  -9.0000000e+01  -9.0000000e+01  -9.0000000e+01  -9.0000000e+01  -9.0000000e+01  -9.0000000e+01  -9.0000000e+01  -3.0000000e+01  -3.0000000e+01  -3.0000000e+01  -3.0000000e+01  -3.0000000e+01  -3.0000000e+01  -3.0000000e+01  -3.0000000e+01   3.0000000e+01   3.0000000e+01   3.0000000e+01   3.0000000e+01   3.0000000e+01   3.0000000e+01   3.0000000e+01   3.0000000e+01   9.0000000e+01   9.0000000e+01   9.0000000e+01   9.0000000e+01   9.0000000e+01   9.0000000e+01   9.0000000e+01   9.0000000e+01   1.5000000e+02   1.5000000e+02   1.5000000e+02   1.5000000e+02   1.5000000e+02   1.5000000e+02   1.5000000e+02   1.5000000e+02   2.1000000e+02   2.1000000e+02   2.1000000e+02   2.1000000e+02   2.1000000e+02   2.1000000e+02   2.1000000e+02   2.1000000e+02   2.7000000e+02   2.7000000e+02   2.7000000e+02   2.7000000e+02   2.7000000e+02   2.7000000e+02   2.7000000e+02   2.7000000e+02   3.3000000e+02   3.3000000e+02   3.3000000e+02   3.3000000e+02   3.3000000e+02   3.3000000e+02   3.3000000e+02   3.3000000e+02   3.9000000e+02   3.9000000e+02   3.9000000e+02   3.9000000e+02   3.9000000e+02   3.9000000e+02   3.9000000e+02   3.9000000e+02   4.5000000e+02   4.5000000e+02   4.5000000e+02   4.5000000e+02   4.5000000e+02   4.5000000e+02   4.5000000e+02   4.5000000e+02   3.9000000e+02   2.7000000e+02   1.5000000e+02   3.0000000e+01   4.5000000e+02   3.3000000e+02   2.1000000e+02   9.0000000e+01   3.9000000e+02   2.7000000e+02   1.5000000e+02   3.0000000e+01   4.5000000e+02   3.3000000e+02   2.1000000e+02   9.0000000e+01   3.9000000e+02   2.7000000e+02   1.5000000e+02   3.0000000e+01   4.5000000e+02   3.3000000e+02   2.1000000e+02   9.0000000e+01   3.9000000e+02   2.7000000e+02   1.5000000e+02   3.0000000e+01   4.5000000e+02   3.3000000e+02   2.1000000e+02   9.0000000e+01   3.9000000e+02   2.7000000e+02   1.5000000e+02   3.0000000e+01   4.5000000e+02   3.3000000e+02   2.1000000e+02   9.0000000e+01   3.9000000e+02   2.7000000e+02   1.5000000e+02   3.0000000e+01   4.5000000e+02   3.3000000e+02   2.1000000e+02   9.0000000e+01   3.9000000e+02   2.7000000e+02   1.5000000e+02   3.0000000e+01   4.5000000e+02   3.3000000e+02   2.1000000e+02   9.0000000e+01   3.9000000e+02   2.7000000e+02   1.5000000e+02   3.0000000e+01   4.5000000e+02   3.3000000e+02   2.1000000e+02   9.0000000e+01   3.9000000e+02   2.7000000e+02   1.5000000e+02   3.0000000e+01   4.5000000e+02   3.3000000e+02   2.1000000e+02   9.0000000e+01   3.9000000e+02   2.7000000e+02   1.5000000e+02   3.0000000e+01   4.5000000e+02   3.3000000e+02   2.1000000e+02   9.0000000e+01   3.9000000e+02   2.7000000e+02   1.5000000e+02   3.0000000e+01   4.5000000e+02   3.3000000e+02   2.1000000e+02   9.0000000e+01   3.9000000e+02   2.7000000e+02   1.5000000e+02   3.0000000e+01   4.5000000e+02   3.3000000e+02   2.1000000e+02   9.0000000e+01   3.9000000e+02   2.7000000e+02   1.5000000e+02   3.0000000e+01   4.5000000e+02   3.3000000e+02   2.1000000e+02   9.0000000e+01   3.9000000e+02   2.7000000e+02   1.5000000e+02   3.0000000e+01   4.5000000e+02   3.3000000e+02   2.1000000e+02   9.0000000e+01   3.9000000e+02   2.7000000e+02   1.5000000e+02   3.0000000e+01   4.5000000e+02   3.3000000e+02   2.1000000e+02   9.0000000e+01   3.9000000e+02   2.7000000e+02   1.5000000e+02   3.0000000e+01   4.5000000e+02   3.3000000e+02   2.1000000e+02   9.0000000e+01   4.5000000e+02   4.5000000e+02   4.5000000e+02   4.5000000e+02   4.5000000e+02   4.5000000e+02   4.5000000e+02   4.5000000e+02   3.9000000e+02   3.9000000e+02   3.9000000e+02   3.9000000e+02   3.9000000e+02   3.9000000e+02   3.9000000e+02   3.9000000e+02   3.3000000e+02   3.3000000e+02   3.3000000e+02   3.3000000e+02   3.3000000e+02   3.3000000e+02   3.3000000e+02   3.3000000e+02   2.7000000e+02   2.7000000e+02   2.7000000e+02   2.7000000e+02   2.7000000e+02   2.7000000e+02   2.7000000e+02   2.7000000e+02   2.1000000e+02   2.1000000e+02   2.1000000e+02   2.1000000e+02   2.1000000e+02   2.1000000e+02   2.1000000e+02   2.1000000e+02   1.5000000e+02   1.5000000e+02   1.5000000e+02   1.5000000e+02   1.5000000e+02   1.5000000e+02   1.5000000e+02   1.5000000e+02   9.0000000e+01   9.0000000e+01   9.0000000e+01   9.0000000e+01   9.0000000e+01   9.0000000e+01   9.0000000e+01   9.0000000e+01   3.0000000e+01   3.0000000e+01   3.0000000e+01   3.0000000e+01   3.0000000e+01   3.0000000e+01   3.0000000e+01   3.0000000e+01  -3.0000000e+01  -3.0000000e+01  -3.0000000e+01  -3.0000000e+01  -3.0000000e+01  -3.0000000e+01  -3.0000000e+01  -3.0000000e+01  -9.0000000e+01  -9.0000000e+01  -9.0000000e+01  -9.0000000e+01  -9.0000000e+01  -9.0000000e+01  -9.0000000e+01  -9.0000000e+01  -1.5000000e+02  -1.5000000e+02  -1.5000000e+02  -1.5000000e+02  -1.5000000e+02  -1.5000000e+02  -1.5000000e+02  -1.5000000e+02  -2.1000000e+02  -2.1000000e+02  -2.1000000e+02  -2.1000000e+02  -2.1000000e+02  -2.1000000e+02  -2.1000000e+02  -2.1000000e+02  -2.7000000e+02  -2.7000000e+02  -2.7000000e+02  -2.7000000e+02  -2.7000000e+02  -2.7000000e+02  -2.7000000e+02  -2.7000000e+02  -3.3000000e+02  -3.3000000e+02  -3.3000000e+02  -3.3000000e+02  -3.3000000e+02  -3.3000000e+02  -3.3000000e+02  -3.3000000e+02  -3.9000000e+02  -3.9000000e+02  -3.9000000e+02  -3.9000000e+02  -3.9000000e+02  -3.9000000e+02  -3.9000000e+02  -3.9000000e+02  -4.5000000e+02  -4.5000000e+02  -4.5000000e+02  -4.5000000e+02  -4.5000000e+02  -4.5000000e+02  -4.5000000e+02  -4.5000000e+02];
end

if params.array_type==519
    temp=[-4.5000000e+02  -3.9000000e+02  -4.2000000e+02  -4.2000000e+02  -4.2000000e+02  -3.9000000e+02  -3.9000000e+02  -2.7000000e+02  -3.3000000e+02  -3.6000000e+02  -3.6000000e+02  -3.6000000e+02  -2.4000000e+02  -1.8000000e+02  -3.3000000e+02  -3.3000000e+02  -3.0000000e+02  -2.7000000e+02  -3.0000000e+02  -3.0000000e+02  -3.0000000e+02  -1.5000000e+02  -2.1000000e+02  -2.7000000e+02  -2.7000000e+02  -2.4000000e+02  -2.1000000e+02  -2.4000000e+02  -2.4000000e+02  -2.4000000e+02  -9.0000000e+01  -1.5000000e+02  -2.1000000e+02  -2.1000000e+02  -1.8000000e+02  -1.5000000e+02  -1.8000000e+02  -1.8000000e+02  -1.8000000e+02  -1.2000000e+02  -9.0000000e+01  -1.5000000e+02  -1.5000000e+02  -1.2000000e+02  -9.0000000e+01  -1.2000000e+02  -1.2000000e+02  -1.2000000e+02  -6.0000000e+01  -3.0000000e+01  -9.0000000e+01  -9.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -3.0000000e+01  -3.0000000e+01  -3.0000000e+01  -3.0000000e+01  -3.0000000e+01   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   3.0000000e+01   3.0000000e+01   3.0000000e+01   3.0000000e+01   3.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   9.0000000e+01   9.0000000e+01   0.0000000e+00   6.0000000e+01   1.2000000e+02   1.2000000e+02   1.2000000e+02   9.0000000e+01   1.2000000e+02   1.5000000e+02   1.5000000e+02   9.0000000e+01   1.2000000e+02   1.8000000e+02   1.8000000e+02   1.8000000e+02   1.5000000e+02   1.8000000e+02   2.1000000e+02   2.1000000e+02   1.5000000e+02   9.0000000e+01   2.4000000e+02   2.4000000e+02   2.4000000e+02   2.1000000e+02   2.4000000e+02   2.7000000e+02   2.7000000e+02   2.1000000e+02   1.5000000e+02   3.0000000e+02   3.0000000e+02   3.0000000e+02   2.7000000e+02   3.0000000e+02   3.3000000e+02   3.3000000e+02   1.8000000e+02   2.4000000e+02   3.6000000e+02   3.6000000e+02   3.6000000e+02   3.3000000e+02   3.9000000e+02   3.9000000e+02   4.2000000e+02   4.2000000e+02   4.2000000e+02   3.9000000e+02   4.5000000e+02   2.7000000e+02   4.8000000e+02   4.5000000e+02   5.1000000e+02   5.4000000e+02   4.8000000e+02   5.7000000e+02   5.1000000e+02   4.5000000e+02   4.8000000e+02   5.4000000e+02   6.0000000e+02   4.2000000e+02   3.6000000e+02   3.3000000e+02   3.9000000e+02   4.5000000e+02   5.1000000e+02   5.7000000e+02   6.3000000e+02   3.0000000e+02   3.6000000e+02   4.2000000e+02   4.8000000e+02   5.4000000e+02   6.0000000e+02   6.6000000e+02   1.2000000e+02   2.1000000e+02   2.7000000e+02   3.3000000e+02   3.9000000e+02   4.5000000e+02   5.1000000e+02   5.7000000e+02   6.3000000e+02   6.9000000e+02   1.8000000e+02   2.4000000e+02   3.0000000e+02   3.6000000e+02   4.2000000e+02   4.8000000e+02   5.4000000e+02   6.0000000e+02   6.6000000e+02   7.2000000e+02   2.1000000e+02   2.7000000e+02   3.3000000e+02   3.9000000e+02   4.5000000e+02   5.1000000e+02   5.7000000e+02   6.3000000e+02   6.9000000e+02   7.5000000e+02   1.5000000e+02   9.0000000e+01   3.0000000e+01   6.0000000e+01   1.8000000e+02   3.0000000e+02   4.2000000e+02   5.4000000e+02   6.6000000e+02   7.8000000e+02   7.2000000e+02   6.0000000e+02   4.8000000e+02   3.6000000e+02   2.4000000e+02   1.2000000e+02   3.0000000e+01   9.0000000e+01   1.5000000e+02   7.5000000e+02   6.9000000e+02   6.3000000e+02   5.7000000e+02   5.1000000e+02   4.5000000e+02   3.9000000e+02   3.3000000e+02   2.7000000e+02   2.1000000e+02   7.2000000e+02   6.6000000e+02   6.0000000e+02   5.4000000e+02   4.8000000e+02   4.2000000e+02   3.6000000e+02   3.0000000e+02   2.4000000e+02   1.8000000e+02   6.9000000e+02   6.3000000e+02   5.7000000e+02   5.1000000e+02   4.5000000e+02   3.9000000e+02   3.3000000e+02   2.7000000e+02   2.1000000e+02   1.2000000e+02   6.6000000e+02   6.0000000e+02   5.4000000e+02   4.8000000e+02   4.2000000e+02   3.6000000e+02   3.0000000e+02   6.3000000e+02   5.7000000e+02   5.1000000e+02   4.5000000e+02   3.9000000e+02   3.3000000e+02   3.6000000e+02   4.2000000e+02   6.0000000e+02   5.4000000e+02   4.8000000e+02   4.5000000e+02   5.1000000e+02   5.7000000e+02   4.8000000e+02   5.4000000e+02   5.1000000e+02   4.8000000e+02   4.5000000e+02   4.5000000e+02   3.9000000e+02   4.2000000e+02   4.2000000e+02   4.2000000e+02   3.9000000e+02   3.9000000e+02   2.7000000e+02   3.3000000e+02   3.6000000e+02   3.6000000e+02   3.6000000e+02   2.4000000e+02   1.8000000e+02   3.3000000e+02   3.3000000e+02   3.0000000e+02   2.7000000e+02   3.0000000e+02   3.0000000e+02   3.0000000e+02   1.5000000e+02   2.1000000e+02   2.7000000e+02   2.7000000e+02   2.4000000e+02   2.1000000e+02   2.4000000e+02   2.4000000e+02   2.4000000e+02   9.0000000e+01   1.5000000e+02   2.1000000e+02   2.1000000e+02   1.8000000e+02   1.5000000e+02   1.8000000e+02   1.8000000e+02   1.8000000e+02   1.2000000e+02   9.0000000e+01   1.5000000e+02   1.5000000e+02   1.2000000e+02   9.0000000e+01   1.2000000e+02   1.2000000e+02   1.2000000e+02   6.0000000e+01   0.0000000e+00   9.0000000e+01   9.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   3.0000000e+01   3.0000000e+01   3.0000000e+01   3.0000000e+01   3.0000000e+01   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00  -3.0000000e+01  -3.0000000e+01  -3.0000000e+01  -3.0000000e+01  -3.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -9.0000000e+01  -9.0000000e+01  -3.0000000e+01  -6.0000000e+01  -1.2000000e+02  -1.2000000e+02  -1.2000000e+02  -9.0000000e+01  -1.2000000e+02  -1.5000000e+02  -1.5000000e+02  -9.0000000e+01  -1.2000000e+02  -1.8000000e+02  -1.8000000e+02  -1.8000000e+02  -1.5000000e+02  -1.8000000e+02  -2.1000000e+02  -2.1000000e+02  -1.5000000e+02  -9.0000000e+01  -2.4000000e+02  -2.4000000e+02  -2.4000000e+02  -2.1000000e+02  -2.4000000e+02  -2.7000000e+02  -2.7000000e+02  -2.1000000e+02  -1.5000000e+02  -3.0000000e+02  -3.0000000e+02  -3.0000000e+02  -2.7000000e+02  -3.0000000e+02  -3.3000000e+02  -3.3000000e+02  -1.8000000e+02  -2.4000000e+02  -3.6000000e+02  -3.6000000e+02  -3.6000000e+02  -3.3000000e+02  -2.7000000e+02  -3.9000000e+02  -3.9000000e+02  -4.2000000e+02  -4.2000000e+02  -4.2000000e+02  -3.9000000e+02  -4.5000000e+02  -4.8000000e+02  -5.1000000e+02  -4.5000000e+02  -5.4000000e+02  -4.8000000e+02  -5.7000000e+02  -5.1000000e+02  -4.5000000e+02  -4.8000000e+02  -5.4000000e+02  -6.0000000e+02  -4.2000000e+02  -3.6000000e+02  -3.3000000e+02  -3.9000000e+02  -4.5000000e+02  -5.1000000e+02  -5.7000000e+02  -6.3000000e+02  -3.0000000e+02  -3.6000000e+02  -4.2000000e+02  -4.8000000e+02  -5.4000000e+02  -6.0000000e+02  -6.6000000e+02  -1.2000000e+02  -2.1000000e+02  -2.7000000e+02  -3.3000000e+02  -3.9000000e+02  -4.5000000e+02  -5.1000000e+02  -5.7000000e+02  -6.3000000e+02  -6.9000000e+02  -1.8000000e+02  -2.4000000e+02  -3.0000000e+02  -3.6000000e+02  -4.2000000e+02  -4.8000000e+02  -5.4000000e+02  -6.0000000e+02  -6.6000000e+02  -7.2000000e+02  -2.1000000e+02  -2.7000000e+02  -3.3000000e+02  -3.9000000e+02  -4.5000000e+02  -5.1000000e+02  -5.7000000e+02  -6.3000000e+02  -6.9000000e+02  -7.5000000e+02  -1.5000000e+02  -9.0000000e+01   0.0000000e+00  -1.2000000e+02  -2.4000000e+02  -3.6000000e+02  -4.8000000e+02  -6.0000000e+02  -7.2000000e+02  -7.8000000e+02  -6.6000000e+02  -5.4000000e+02  -4.2000000e+02  -3.0000000e+02  -1.8000000e+02  -6.0000000e+01  -9.0000000e+01  -1.5000000e+02  -7.5000000e+02  -6.9000000e+02  -6.3000000e+02  -5.7000000e+02  -5.1000000e+02  -4.5000000e+02  -3.9000000e+02  -3.3000000e+02  -2.7000000e+02  -2.1000000e+02  -7.2000000e+02  -6.6000000e+02  -6.0000000e+02  -5.4000000e+02  -4.8000000e+02  -4.2000000e+02  -3.6000000e+02  -3.0000000e+02  -2.4000000e+02  -1.8000000e+02  -6.9000000e+02  -6.3000000e+02  -5.7000000e+02  -5.1000000e+02  -4.5000000e+02  -3.9000000e+02  -3.3000000e+02  -2.7000000e+02  -2.1000000e+02  -1.2000000e+02  -6.6000000e+02  -6.0000000e+02  -5.4000000e+02  -4.8000000e+02  -4.2000000e+02  -3.6000000e+02  -3.0000000e+02  -6.3000000e+02  -5.7000000e+02  -5.1000000e+02  -4.5000000e+02  -3.9000000e+02  -3.3000000e+02  -3.6000000e+02  -4.2000000e+02  -6.0000000e+02  -5.4000000e+02  -4.8000000e+02  -4.5000000e+02  -5.1000000e+02  -5.7000000e+02  -4.8000000e+02  -5.4000000e+02  -4.5000000e+02  -5.1000000e+02  -4.8000000e+02;...
        -6.6000000e+02  -4.2000000e+02  -4.8000000e+02  -7.2000000e+02  -6.0000000e+02  -5.4000000e+02  -6.6000000e+02  -3.0000000e+02  -4.2000000e+02  -4.8000000e+02  -7.2000000e+02  -6.0000000e+02  -2.4000000e+02  -2.4000000e+02  -5.4000000e+02  -6.6000000e+02  -3.6000000e+02  -4.2000000e+02  -4.8000000e+02  -7.2000000e+02  -6.0000000e+02  -1.8000000e+02  -3.0000000e+02  -5.4000000e+02  -6.6000000e+02  -3.6000000e+02  -4.2000000e+02  -4.8000000e+02  -7.2000000e+02  -6.0000000e+02  -1.8000000e+02  -3.0000000e+02  -5.4000000e+02  -6.6000000e+02  -3.6000000e+02  -4.2000000e+02  -4.8000000e+02  -7.2000000e+02  -6.0000000e+02  -2.4000000e+02  -3.0000000e+02  -5.4000000e+02  -6.6000000e+02  -3.6000000e+02  -4.2000000e+02  -4.8000000e+02  -7.2000000e+02  -6.0000000e+02  -1.2000000e+02  -6.0000000e+01  -5.4000000e+02  -6.6000000e+02  -2.4000000e+02  -3.6000000e+02  -4.8000000e+02  -7.2000000e+02  -6.0000000e+02  -1.8000000e+02  -3.0000000e+02  -5.4000000e+02  -6.6000000e+02  -4.2000000e+02  -2.4000000e+02  -4.8000000e+02  -7.2000000e+02  -6.0000000e+02  -3.6000000e+02  -4.2000000e+02  -6.6000000e+02  -5.4000000e+02  -3.0000000e+02  -1.8000000e+02  -6.0000000e+02  -7.2000000e+02  -4.8000000e+02  -3.6000000e+02  -2.4000000e+02  -6.6000000e+02  -5.4000000e+02  -1.2000000e+02  -1.2000000e+02  -6.0000000e+02  -7.2000000e+02  -4.8000000e+02  -4.2000000e+02  -3.6000000e+02  -6.6000000e+02  -5.4000000e+02  -3.0000000e+02  -2.4000000e+02  -6.0000000e+02  -7.2000000e+02  -4.8000000e+02  -4.2000000e+02  -3.6000000e+02  -6.6000000e+02  -5.4000000e+02  -3.0000000e+02  -1.8000000e+02  -6.0000000e+02  -7.2000000e+02  -4.8000000e+02  -4.2000000e+02  -3.6000000e+02  -6.6000000e+02  -5.4000000e+02  -3.0000000e+02  -1.8000000e+02  -6.0000000e+02  -7.2000000e+02  -4.8000000e+02  -4.2000000e+02  -3.6000000e+02  -6.6000000e+02  -5.4000000e+02  -2.4000000e+02  -2.4000000e+02  -6.0000000e+02  -7.2000000e+02  -4.8000000e+02  -4.2000000e+02  -6.6000000e+02  -5.4000000e+02  -6.0000000e+02  -7.2000000e+02  -4.8000000e+02  -4.2000000e+02  -6.6000000e+02  -3.0000000e+02  -6.0000000e+02  -5.4000000e+02  -5.4000000e+02  -4.8000000e+02  -4.8000000e+02  -4.2000000e+02  -4.2000000e+02  -4.2000000e+02  -3.6000000e+02  -3.6000000e+02  -3.6000000e+02  -3.6000000e+02  -3.6000000e+02  -3.0000000e+02  -3.0000000e+02  -3.0000000e+02  -3.0000000e+02  -3.0000000e+02  -3.0000000e+02  -2.4000000e+02  -2.4000000e+02  -2.4000000e+02  -2.4000000e+02  -2.4000000e+02  -2.4000000e+02  -2.4000000e+02  -1.2000000e+02  -1.8000000e+02  -1.8000000e+02  -1.8000000e+02  -1.8000000e+02  -1.8000000e+02  -1.8000000e+02  -1.8000000e+02  -1.8000000e+02  -1.8000000e+02  -1.2000000e+02  -1.2000000e+02  -1.2000000e+02  -1.2000000e+02  -1.2000000e+02  -1.2000000e+02  -1.2000000e+02  -1.2000000e+02  -1.2000000e+02  -1.2000000e+02  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   6.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   1.2000000e+02   1.2000000e+02   1.2000000e+02   1.2000000e+02   1.2000000e+02   1.2000000e+02   1.2000000e+02   1.2000000e+02   1.2000000e+02   1.2000000e+02   1.8000000e+02   1.8000000e+02   1.8000000e+02   1.8000000e+02   1.8000000e+02   1.8000000e+02   1.8000000e+02   1.8000000e+02   1.8000000e+02   1.2000000e+02   2.4000000e+02   2.4000000e+02   2.4000000e+02   2.4000000e+02   2.4000000e+02   2.4000000e+02   2.4000000e+02   3.0000000e+02   3.0000000e+02   3.0000000e+02   3.0000000e+02   3.0000000e+02   3.0000000e+02   3.6000000e+02   3.6000000e+02   3.6000000e+02   3.6000000e+02   3.6000000e+02   4.2000000e+02   4.2000000e+02   4.2000000e+02   4.8000000e+02   4.8000000e+02   5.4000000e+02   6.0000000e+02   5.4000000e+02   6.6000000e+02   4.2000000e+02   4.8000000e+02   7.2000000e+02   6.0000000e+02   5.4000000e+02   6.6000000e+02   3.0000000e+02   4.2000000e+02   4.8000000e+02   7.2000000e+02   6.0000000e+02   2.4000000e+02   2.4000000e+02   5.4000000e+02   6.6000000e+02   3.6000000e+02   4.2000000e+02   4.8000000e+02   7.2000000e+02   6.0000000e+02   1.8000000e+02   3.0000000e+02   5.4000000e+02   6.6000000e+02   3.6000000e+02   4.2000000e+02   4.8000000e+02   7.2000000e+02   6.0000000e+02   1.8000000e+02   3.0000000e+02   5.4000000e+02   6.6000000e+02   3.6000000e+02   4.2000000e+02   4.8000000e+02   7.2000000e+02   6.0000000e+02   2.4000000e+02   3.0000000e+02   5.4000000e+02   6.6000000e+02   3.6000000e+02   4.2000000e+02   4.8000000e+02   7.2000000e+02   6.0000000e+02   1.2000000e+02   1.2000000e+02   5.4000000e+02   6.6000000e+02   2.4000000e+02   3.6000000e+02   4.8000000e+02   7.2000000e+02   6.0000000e+02   1.8000000e+02   3.0000000e+02   5.4000000e+02   6.6000000e+02   4.2000000e+02   3.6000000e+02   6.0000000e+02   7.2000000e+02   4.8000000e+02   2.4000000e+02   4.2000000e+02   6.6000000e+02   5.4000000e+02   3.0000000e+02   1.8000000e+02   6.0000000e+02   7.2000000e+02   4.8000000e+02   3.6000000e+02   2.4000000e+02   6.6000000e+02   5.4000000e+02   6.0000000e+01   1.2000000e+02   6.0000000e+02   7.2000000e+02   4.8000000e+02   4.2000000e+02   3.6000000e+02   6.6000000e+02   5.4000000e+02   3.0000000e+02   2.4000000e+02   6.0000000e+02   7.2000000e+02   4.8000000e+02   4.2000000e+02   3.6000000e+02   6.6000000e+02   5.4000000e+02   3.0000000e+02   1.8000000e+02   6.0000000e+02   7.2000000e+02   4.8000000e+02   4.2000000e+02   3.6000000e+02   6.6000000e+02   5.4000000e+02   3.0000000e+02   1.8000000e+02   6.0000000e+02   7.2000000e+02   4.8000000e+02   4.2000000e+02   3.6000000e+02   6.6000000e+02   5.4000000e+02   2.4000000e+02   2.4000000e+02   6.0000000e+02   7.2000000e+02   4.8000000e+02   4.2000000e+02   3.0000000e+02   6.6000000e+02   5.4000000e+02   6.0000000e+02   7.2000000e+02   4.8000000e+02   4.2000000e+02   6.6000000e+02   6.0000000e+02   5.4000000e+02   5.4000000e+02   4.8000000e+02   4.8000000e+02   4.2000000e+02   4.2000000e+02   4.2000000e+02   3.6000000e+02   3.6000000e+02   3.6000000e+02   3.6000000e+02   3.6000000e+02   3.0000000e+02   3.0000000e+02   3.0000000e+02   3.0000000e+02   3.0000000e+02   3.0000000e+02   2.4000000e+02   2.4000000e+02   2.4000000e+02   2.4000000e+02   2.4000000e+02   2.4000000e+02   2.4000000e+02   1.2000000e+02   1.8000000e+02   1.8000000e+02   1.8000000e+02   1.8000000e+02   1.8000000e+02   1.8000000e+02   1.8000000e+02   1.8000000e+02   1.8000000e+02   1.2000000e+02   1.2000000e+02   1.2000000e+02   1.2000000e+02   1.2000000e+02   1.2000000e+02   1.2000000e+02   1.2000000e+02   1.2000000e+02   1.2000000e+02   6.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   6.0000000e+01   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00   0.0000000e+00  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -6.0000000e+01  -1.2000000e+02  -1.2000000e+02  -1.2000000e+02  -1.2000000e+02  -1.2000000e+02  -1.2000000e+02  -1.2000000e+02  -1.2000000e+02  -1.2000000e+02  -1.2000000e+02  -1.8000000e+02  -1.8000000e+02  -1.8000000e+02  -1.8000000e+02  -1.8000000e+02  -1.8000000e+02  -1.8000000e+02  -1.8000000e+02  -1.8000000e+02  -1.2000000e+02  -2.4000000e+02  -2.4000000e+02  -2.4000000e+02  -2.4000000e+02  -2.4000000e+02  -2.4000000e+02  -2.4000000e+02  -3.0000000e+02  -3.0000000e+02  -3.0000000e+02  -3.0000000e+02  -3.0000000e+02  -3.0000000e+02  -3.6000000e+02  -3.6000000e+02  -3.6000000e+02  -3.6000000e+02  -3.6000000e+02  -4.2000000e+02  -4.2000000e+02  -4.2000000e+02  -4.8000000e+02  -4.8000000e+02  -5.4000000e+02  -5.4000000e+02  -6.0000000e+02];
end

datarun.position=temp';