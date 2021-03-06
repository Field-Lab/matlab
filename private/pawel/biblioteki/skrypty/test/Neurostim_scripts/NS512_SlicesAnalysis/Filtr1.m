function Hd = Filtr1
%FILTR1 Returns a discrete-time filter object.

%
% M-File generated by MATLAB(R) 7.9 and the Signal Processing Toolbox 6.12.
%
% Generated on: 26-Jun-2014 14:55:24
%

Fstop = 300;    % Stopband Frequency
Fpass = 500;    % Passband Frequency
Astop = 50;     % Stopband Attenuation (dB)
Apass = 0.5;    % Passband Ripple (dB)
Fs    = 20000;  % Sampling Frequency

h = fdesign.highpass('fst,fp,ast,ap', Fstop, Fpass, Astop, Apass, Fs);

Hd = design(h, 'equiripple', ...
    'MinOrder', 'any');



% [EOF]
