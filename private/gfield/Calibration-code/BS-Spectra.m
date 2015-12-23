clear all

load red
load green
load blue

RGBMatrix = [red(1:200) green(1:200) blue(1:200)];

save RGBMatrix RGBMatrix

save rgb.2 RGBMatrix -ascii


%%%%%
load RigALargeBeamSplitter
load RigASpectra

Ri

LargeBeamSplitterSpectrum = RigALargeBeamSplitter ./ RigASpectra;

figure(2)
plot(LargeBeamSplitterSpectrum)
hold on
plot(0.5*ones(1,200), 'r')
hold off


NoiseSamplesIndex = find(LargeBeamSplitterSpectrum > 0.55);

SmoothedBeamSplitterSpectrum = LargeBeamSplitterSpectrum;
SmoothedBeamSplitterSpectrum(NoiseSamplesIndex) = 0.55;

figure(1)
plot(SmoothedBeamSplitterSpectrum)


load A-rig
load BSA

BS = BSA ./ A_rig;
plot(BS)
hold on

%%%%%%
load RigASpec
load LargeBSSpec
load SmallBSSpec

BSFilter = SmallBSSpec ./ RigASpec
figure(1)
plot(BSFilter)

Threshold = 0.485
ThresholdIndex = find(BSFilter > Threshold);

ThresholdedFilter = BSFilter;
ThresholdedFilter(ThresholdIndex(1):length(BSFilter)) = Threshold

figure(2)
plot(ThresholdedFilter)
axis([1 200 0.2 0.6])

save LargeBS ThresholdedFilter




