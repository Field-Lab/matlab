function Hd = filtrhigh
%FILTRHIGH Returns a discrete-time filter object.

%
% M-File generated by MATLAB(R) 7.9 and the Signal Processing Toolbox 6.12.
%
% Generated on: 06-May-2013 14:40:41
%

% Equiripple Highpass filter designed using the FIRPM function.

% All frequency values are in Hz.
Fs = 20000;  % Sampling Frequency

N     = 201;  % Order
Fstop = 20;   % Stopband Frequency
Fpass = 30;   % Passband Frequency
Wstop = 1;    % Stopband Weight
Wpass = 1;    % Passband Weight
dens  = 20;   % Density Factor

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, [0 Fstop Fpass Fs/2]/(Fs/2), [0 0 1 1], [Wstop Wpass], ...
           {dens});
Hd = dfilt.dffir(b);

% [EOF]
