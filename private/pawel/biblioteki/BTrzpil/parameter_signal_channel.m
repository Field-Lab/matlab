%parametry impulsu
%amplitude kolejne amplitudy impulsu
%time czas trwania poszczeg�lnych amplitud wyrazone w liczbach calkowitych
%- wielokrotnosc od�wie�ania
%channel - kanal
%start - ktora wielokrotno�� czasu od�wiezania jest poczatkiem impulsu
function [Amplituda, Range, T_amp, T_time, channel, start]=parameter_signal_channel(A, range, amp, t, ch, str)
T_amp=[amp];
T_time=[t];
channel=ch;
start=str;
Range=range;
Amplituda=A;
