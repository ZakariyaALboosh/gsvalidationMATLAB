function [snr_meas_dB, track, time_offset_s] = extractMeasuredSNR(wf, geo, opts)
% EXTRACTMEASUREDSNR  Extract per-row SNR from a waterfall along the Doppler track.
%
%   For each waterfall row: a signal window of width sig_bw_Hz is centred
%   on the expected signal track; S = mean of the top-3 bins in that
%   window (dB), N = median of the noise-guard bins (dB), and
%   SNR_meas = S - N. Because both S and N are read from the SAME
%   uncalibrated (or std-scaled) waterfall, any unknown gain/scaling
%   constant cancels in the difference — this is why the pipeline
%   validates SNR rather than absolute power.
%
%   TRACK MODEL: track(t) = doppler_factor * doppler(t) + freq_offset_Hz.
%   SatNOGS waterfalls (network .h5 artifacts and the client .dat they come
%   from) are generated AFTER the client's Doppler correction (verified
%   empirically: a fixed terrestrial carrier in a real artifact follows
%   -doppler(t) to within one bin — the exact mirror of the correction).
%   So for SatNOGS waterfalls use doppler_factor = 0: the satellite sits
%   near a CONSTANT offset (transmitter offset + oscillator drift), and
%   fixed local carriers appear as inverted S-curves. doppler_factor = 1
%   is for genuinely uncorrected spectra.
%
% Inputs:
%   wf   - waterfall struct from parseWaterfall / generateSyntheticWaterfall:
%            .t_s (Nx1 s), .f_Hz (1xM Hz offsets), .P (NxM dB-like),
%            .start_time (datetime UTC or NaT)
%   geo  - geometry struct: .time (datetime UTC), .doppler_Hz (Hz)
%   opts - optional struct, any subset of:
%            .sig_bw_Hz     signal window width (Hz)          [default 5000]
%            .guard_lo_Hz   inner edge of noise guard band, offset from
%                           track centre (Hz)                 [default 20e3]
%            .guard_hi_Hz   outer edge of noise guard band (Hz) [default 50e3]
%            .time_offset_s manual nudge added to the waterfall->geometry
%                           time alignment (s)                [default 0]
%            .doppler_factor multiplier on the predicted Doppler in the
%                           track model: 1 = uncorrected spectrum,
%                           0 = Doppler-corrected (SatNOGS)    [default 1]
%            .freq_offset_Hz constant track offset, e.g. transmitter
%                           off-frequency (Hz)                [default 0]
%
% Outputs:
%   snr_meas_dB   - Nx1 measured SNR per waterfall row (dB). NaN where the
%                   row falls outside the geometry time span.
%   track         - Nx1 frequency-bin index of the signal-track centre per
%                   row (NaN where undefined); for diagnostic overlay on
%                   the waterfall image.
%   time_offset_s - scalar: the effective offset (s) that was
%                   ADDED to wf.t_s to place rows on the geo.time axis
%                   (start-time alignment + opts.time_offset_s).

if nargin < 3, opts = struct(); end
if ~isfield(opts, 'sig_bw_Hz'),     opts.sig_bw_Hz     = 5000; end
if ~isfield(opts, 'guard_lo_Hz'),   opts.guard_lo_Hz   = 20e3; end
if ~isfield(opts, 'guard_hi_Hz'),   opts.guard_hi_Hz   = 50e3; end
if ~isfield(opts, 'time_offset_s'),  opts.time_offset_s  = 0;  end
if ~isfield(opts, 'doppler_factor'), opts.doppler_factor = 1;  end
if ~isfield(opts, 'freq_offset_Hz'), opts.freq_offset_Hz = 0;  end

% ---- time alignment: waterfall rows -> geometry timeline --------------
% If the waterfall carries an absolute start_time, anchor to it; otherwise
% assume the waterfall starts at geo.time(1). opts.time_offset_s nudges
% either case.
geo_t_s = seconds(geo.time(:) - geo.time(1));
if isfield(wf, 'start_time') && isdatetime(wf.start_time) && ~isnat(wf.start_time)
    time_offset_s = seconds(wf.start_time - geo.time(1)) + opts.time_offset_s;
else
    time_offset_s = opts.time_offset_s;
end
t_wf_on_geo = wf.t_s(:) + time_offset_s;

dopp_Hz = interp1(geo_t_s, geo.doppler_Hz(:), t_wf_on_geo, 'linear', NaN);
% expected signal position per row (see TRACK MODEL in the header);
% NaN doppler (row outside the geometry span) stays NaN even for factor 0,
% keeping the time masking identical in both modes
track_Hz = opts.doppler_factor * dopp_Hz + opts.freq_offset_Hz;

% ---- per-row extraction ------------------------------------------------
nrows = numel(wf.t_s);
f_Hz  = wf.f_Hz(:).';                     % force 1xM
fmax  = max(abs(f_Hz));                   % passband half-width (Hz)
snr_meas_dB = nan(nrows, 1);
track       = nan(nrows, 1);
guard_warned = false;

for i = 1:nrows
    if ~isfinite(track_Hz(i)), continue; end  % row outside geometry span

    df = f_Hz - track_Hz(i);                  % bin offsets from track centre

    % signal window: +/- sig_bw/2 around the Doppler-shifted carrier
    sig_idx = find(abs(df) <= opts.sig_bw_Hz / 2);
    if isempty(sig_idx), continue; end        % track has left the passband

    % noise guard: annulus guard_lo..guard_hi either side of the track.
    % IMPORTANT: restricted to the inner 80% of the passband, because real
    % SDR waterfalls roll off ~10-30 dB at the band edges (anti-alias
    % filter) - guard bins out there would fake a low noise floor and
    % inflate the SNR (seen on a real SatNOGS artifact: +16 dB bias).
    usable = abs(f_Hz) <= 0.8 * fmax;
    noise_idx = find(abs(df) >= opts.guard_lo_Hz & ...
                     abs(df) <= opts.guard_hi_Hz & usable);
    if numel(noise_idx) < 64        % sparse guard -> noisy median; fall back
        % Guard band unusable (narrow waterfall or guard beyond the flat
        % passband): fall back to the central 60% of the passband, clear
        % of the signal window. The median is robust to the narrow signal
        % and stray carriers in that region.
        noise_idx = find(abs(f_Hz) <= 0.6 * fmax & ...
                         abs(df) > opts.sig_bw_Hz * 1.5);
        if ~guard_warned
            warning('extractMeasuredSNR:guardFallback', ...
                ['Noise guard %g-%g kHz falls outside the flat passband ', ...
                 '(|f| <= %g kHz); using the central-passband median instead.'], ...
                opts.guard_lo_Hz/1e3, opts.guard_hi_Hz/1e3, 0.8*fmax/1e3);
            guard_warned = true;
        end
        if numel(noise_idx) < 8, continue; end   % hopeless: skip the row
    end

    row = wf.P(i, :);
    sig_sorted = sort(row(sig_idx), 'descend');
    S = mean(sig_sorted(1:min(3, numel(sig_sorted))));   % top-3 bins (dB)
    N = median(row(noise_idx));                          % noise floor (dB)
    snr_meas_dB(i) = S - N;

    [~, ic] = min(abs(df));
    track(i) = ic;
end
end
