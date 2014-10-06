function Artifact=NS_PredictArtifact(Amplitude,Art1,Amp1,Art2,Amp2);
%Predicts artifact shape for pulse of given Amplitude based on the artifact
%shape for two different amplitudes. It is assumed that the pulse shape is
%identical for the pulse for which we are prediting the artifact and for
%both pulses for which the artifact is given as input (Art1,Art2).
%This model assumes that the artifact for any amplitude is a combination of
%some waveform which is constant (does not depend on amplitude) and other
%waveform that scales with amplitude in linear fashion (as for 2008-03-15).
%Input:
%Amplitude - amplitude of the pulse for which we need the predicted
%artifact.
%Art1 - artifact for pulse of amplitude Amp1;
%Art2 - artifact for pulse of amplitude Amp2;

Artifact=Art2+(Art2-Art1)*(Amplitude-Amp2)/(Amp2-Amp1); %simple, isn't it