%parametry impulsu
%amplitude kolejne amplitudy impulsu
%time czas trwania poszczególnych amplitud wyrazone w liczbach calkowitych
%- wielokrotnosc odœwie¿ania
%channel - kanal
%start - ktora wielokrotnoœæ czasu odœwiezania jest poczatkiem impulsu
function [Amplituda, Range, T_amp, T_time, channel, start]=parameter_signal_channel(A, range, amp, t, ch, str)
T_amp=[amp];
T_time=[t];
channel=ch;
start=str;
Range=range;
Amplituda=A;
