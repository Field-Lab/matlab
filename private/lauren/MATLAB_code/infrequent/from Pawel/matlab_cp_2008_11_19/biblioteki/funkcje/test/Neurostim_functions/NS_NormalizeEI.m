function NormEI=NS_NormalizeEI(EI);

[Shift,MaxCorr]=NS_FindShiftBetweenEIs(EI,EI);
NormEI=EI/sqrt(MaxCorr);
