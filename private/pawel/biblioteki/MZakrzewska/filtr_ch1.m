function Hd = filtr_ch1
%FILTR_CH1 Returns a discrete-time filter object.

%
% M-File generated by MATLAB(R) 7.9 and the Signal Processing Toolbox 6.12.
%
% Generated on: 08-Aug-2012 11:11:05
%

% Chebyshev Type II Bandstop filter designed using FDESIGN.BANDSTOP.

% All frequency values are in Hz.
Fs = 20000;  % Sampling Frequency

N      = 200;  % Order
Fstop1 = 59;   % First Stopband Frequency
Fstop2 = 61;   % Second Stopband Frequency
Astop  = 40;   % Stopband Attenuation (dB)

% Construct an FDESIGN object and call its CHEBY2 method.
h  = fdesign.bandstop('N,Fst1,Fst2,Ast', N, Fstop1, Fstop2, Astop, Fs);
Hd = design(h, 'cheby2');

% [EOF]
