function powers = extract_fourier_components(trace, deltat, componentfs, sigmafs, varargin)
% Just extracts the power of the requested frequency components
% Defaults to subtract the median power as an attempted baseline correction

opts = inputParser();
opts.addParamValue('subtractmedian', true);
opts.addParamValue('NFFT', []);
opts.parse(varargin{:});
opts = opts.Results;

% Calculate frequencies we expect to get
L = length(trace);
NFFT = opts.NFFT;
if isempty(NFFT), NFFT = 2^nextpow2(L); end
F = 1/deltat;
f = F/2 * linspace(0, 1, NFFT/2+1);

% Find windows around desired frequency components
for i = 1:length(componentfs)
    if length(sigmafs) == 1, sigma = sigmafs;
    else                     sigma = sigmafs(i); end

    fxs{i} = find(abs(f-componentfs(i)) < sigma);
end

% Calculate power spectrum
tracef = fft(trace, NFFT) / L;
pwr = 2 * abs(tracef(1:NFFT/2+1));

% Use median as baseline to subtract?
baseline = 0;
if opts.subtractmedian, baseline = median(pwr); end

% Extract components from windows calculated above
for i = 1:length(fxs)
    powers(i) = sum(pwr(fxs{i})) - baseline*length(fxs{i});
end