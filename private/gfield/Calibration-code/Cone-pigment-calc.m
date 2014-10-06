MolarConcentrationOfPigment = 0.003; %Mols.
AvagadrosNumber = 6e23;
VolumeOfCone = 30; % micron^3 (Schnapf et al 1990)

% includes conversion from L to mL.
MoleculesPer_mL = (MolarConcentrationOfPigment * AvagadrosNumber * 0.001); 
% convert from mL to �m^3 -- 1e4 coverts from cm to �m, for 1 dimension.
MoleculesPerMicronCubed = MoleculesPer_mL / (1e4 * 1e4 * 1e4);
MoleculesPerConeOS = MoleculesPerMicronCubed * VolumeOfCone

% 5.5e7 = molecules per cone
