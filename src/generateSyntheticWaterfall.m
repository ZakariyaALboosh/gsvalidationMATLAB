function wf = generateSyntheticWaterfall(geo, snr_true_dB, params)
% GENERATESYNTHETICWATERFALL  Synthesise a waterfall struct with a known SNR.
%
%   Test harness for the measurement side of the pipeline. Builds a
%   waterfall in the exact field layout returned by parseWaterfall:
%   a Gaussian noise floor (in dB units, mimicking an averaged,
%   uncalibrated spectrum) plus a Doppler-shifted signal ridge whose
%   level is set so that the extraction definition SNR = S - N
%   (top-of-ridge dB minus noise-median dB) recovers snr_true_dB exactly,
%   up to the injected noise jitter.
%
% Inputs:
%   geo         - struct from computePassGeometry (or synthetic equivalent):
%                   .time       datetime vector, UTC (Nx1)
%                   .el_deg     elevation (deg, Nx1)
%                   .doppler_Hz Doppler shift at downlink freq (Hz, Nx1)
%   snr_true_dB - true SNR profile on geo.time (dB, Nx1). Rows where the
%                 profile is NaN (or elevation < 0 deg) get no signal ridge.
%   params      - optional struct, any subset of:
%                   .samp_rate_Hz    waterfall span (Hz)      [default 48e3]
%                   .nchan           frequency bins           [default 1024]
%                   .dt_s            row spacing (s)          [default 1]
%                   .center_Hz       RF centre frequency (Hz) [default 435e6]
%                   .noise_floor_dB  noise floor level (dB, arbitrary ref)
%                                                             [default -100]
%                   .noise_sigma_dB  per-bin noise std (dB)   [default 0.7]
%                   .ridge_bins      width of signal ridge (bins, odd)
%                                                             [default 3]
%
% Outputs:
%   wf - waterfall struct (same layout as parseWaterfall):
%          .t_s        Nrows x 1, seconds since first row
%          .f_Hz       1 x nchan, bin offsets from centre (Hz)
%          .center_Hz  scalar (Hz)
%          .P          Nrows x nchan power (dB, arbitrary reference)
%          .start_time datetime UTC of first row (= geo.time(1), so that
%                      time alignment in extractMeasuredSNR is exact)
%          .meta       empty struct (synthetic data has no artifact metadata)

% ---- defaults --------------------------------------------------------
if nargin < 3, params = struct(); end
p = params;
if ~isfield(p, 'samp_rate_Hz'),   p.samp_rate_Hz   = 48e3;  end
if ~isfield(p, 'nchan'),          p.nchan          = 1024;  end
if ~isfield(p, 'dt_s'),           p.dt_s           = 1;     end
if ~isfield(p, 'center_Hz'),      p.center_Hz      = 435e6; end
if ~isfield(p, 'noise_floor_dB'), p.noise_floor_dB = -100;  end
if ~isfield(p, 'noise_sigma_dB'), p.noise_sigma_dB = 0.7;   end
if ~isfield(p, 'ridge_bins'),     p.ridge_bins     = 3;     end

geo_t_s = seconds(geo.time(:) - geo.time(1));   % relative pass time (s)
snr_true_dB = snr_true_dB(:);
if numel(snr_true_dB) ~= numel(geo_t_s)
    error('generateSyntheticWaterfall:badInput', ...
        'snr_true_dB (%d) must have one value per geo.time sample (%d).', ...
        numel(snr_true_dB), numel(geo_t_s));
end

% ---- axes (same conventions as parseWaterfallDAT) --------------------
t_s  = (0:p.dt_s:geo_t_s(end)).';                             % Nrows x 1
f_Hz = (-0.5 : 1/p.nchan : 0.5 - 1/p.nchan) * p.samp_rate_Hz; % 1 x nchan
nrows = numel(t_s);

% ---- interpolate truth onto the waterfall row times ------------------
dopp_Hz = interp1(geo_t_s, geo.doppler_Hz(:), t_s, 'linear', NaN);
el_deg  = interp1(geo_t_s, geo.el_deg(:),     t_s, 'linear', NaN);
snr_dB  = interp1(geo_t_s, snr_true_dB,       t_s, 'linear', NaN);

% ---- noise floor ------------------------------------------------------
% Gaussian in the dB domain: real waterfall rows are averages of many FFTs
% (nfft_per_row), which shrinks the exponential-power skew toward Gaussian.
P = p.noise_floor_dB + p.noise_sigma_dB * randn(nrows, p.nchan);

% ---- Doppler-shifted signal ridge -------------------------------------
% Ridge bins are set to floor + SNR: since extractMeasuredSNR defines
% SNR_meas = (top-of-ridge dB) - (noise median dB), this makes the
% target profile recoverable exactly by that definition.
half = floor(p.ridge_bins / 2);
for i = 1:nrows
    if ~isfinite(snr_dB(i)) || ~isfinite(dopp_Hz(i)) || el_deg(i) < 0
        continue;                        % satellite below horizon: noise only
    end
    [~, ic] = min(abs(f_Hz - dopp_Hz(i)));
    idx = max(1, ic - half) : min(p.nchan, ic + half);
    P(i, idx) = p.noise_floor_dB + snr_dB(i) ...
                + 0.1 * p.noise_sigma_dB * randn(1, numel(idx));
end

% ---- assemble unified waterfall struct --------------------------------
wf = struct();
wf.t_s        = t_s;
wf.f_Hz       = f_Hz;
wf.center_Hz  = p.center_Hz;
wf.P          = P;
wf.start_time = geo.time(1);
wf.meta       = struct();
end
