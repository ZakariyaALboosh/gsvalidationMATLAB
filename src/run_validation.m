% RUN_VALIDATION  Top-level script for the SatNOGS link-budget validation.
%
%   The ONLY file with hardcoded paths / run configuration. Edit the config
%   block, press F5. Modes:
%     "selftest" - synthetic acceptance test of the SNR extraction and
%                  comparison chain (no toolboxes, no data files needed):
%                  known SNR profile -> generateSyntheticWaterfall ->
%                  extractMeasuredSNR -> compareAndPlot. PASSES when the
%                  extracted SNR matches the input within 0.5 dB (max abs
%                  error above the elevation mask) and correlation > 0.99.
%     "hdf5"     - full pipeline on a SatNOGS network .h5 artifact
%                  (parse -> config -> geometry -> predict -> extract ->
%                  compare).                       [not yet built: step 3/5]
%     "dat"      - full pipeline on a raw client .dat waterfall.
%                                                  [not yet built: step 4/5]

% ======================= CONFIG BLOCK ==================================
MODE = "selftest";            % "selftest" | "hdf5" | "dat"

paths.repo    = fileparts(mfilename('fullpath'));      % .../src
paths.data    = fullfile(paths.repo, '..', 'data');
paths.output  = fullfile(paths.repo, '..', 'output');
paths.config  = fullfile(paths.repo, '..', 'config', 'station_uhf.m');
paths.wf_file = fullfile(paths.data, 'observation.h5'); % .h5 or .dat to run

obs_start     = datetime(2026, 1, 1, 12, 0, 0, 'TimeZone', 'UTC'); % pass window
obs_stop      = obs_start + minutes(12);
time_offset_s = 0;            % manual waterfall<->geometry alignment nudge (s)
% =======================================================================

addpath(paths.repo);
rng(42);                      % reproducible selftest noise

switch MODE
    case "selftest"
        fprintf('=== SELFTEST: synthetic waterfall acceptance test ===\n');

        % --- synthetic pass geometry (no toolbox needed) ---------------
        geo = makeSyntheticGeo(obs_start, 600);

        % --- known "true" SNR profile: rises and falls with elevation --
        snr_true_dB = 4 + 21 * sind(max(geo.el_deg, 0));   % ~4..25 dB
        snr_true_dB(geo.el_deg < 0) = NaN;                 % below horizon

        % --- synthesise -> extract -> compare ---------------------------
        wf = generateSyntheticWaterfall(geo, snr_true_dB);
        [snr_meas_dB, track, toff] = extractMeasuredSNR(wf, geo, ...
            struct('time_offset_s', time_offset_s));
        results = compareAndPlot(geo, snr_true_dB, wf, snr_meas_dB, ...
            struct('out_dir', paths.output, 'label', 'selftest', ...
                   'time_offset_s', toff, 'track', track));

        % --- acceptance criteria ----------------------------------------
        max_abs_err_dB = max(abs(results.snr_meas_dB(results.mask) ...
                                 - results.snr_pred_i_dB(results.mask)));
        fprintf('Selftest: max |meas - true| = %.3f dB (limit 0.50), ', ...
            max_abs_err_dB);
        fprintf('correlation = %.4f (limit 0.99)\n', results.corr);
        assert(max_abs_err_dB <= 0.5, ...
            'SELFTEST FAILED: extracted SNR deviates %.3f dB (> 0.5 dB).', ...
            max_abs_err_dB);
        assert(results.corr > 0.99, ...
            'SELFTEST FAILED: correlation %.4f <= 0.99.', results.corr);
        fprintf('SELFTEST PASSED. Figures and results in %s\n', paths.output);

    case "hdf5"
        error('run_validation:notBuilt', ...
            'hdf5 mode arrives with build steps 3/5 (parseWaterfallHDF5 + integration).');

    case "dat"
        error('run_validation:notBuilt', ...
            'dat mode arrives with build step 4 (parseWaterfallDAT).');

    otherwise
        error('run_validation:badMode', 'Unknown MODE "%s".', MODE);
end

% ======================= LOCAL FUNCTIONS ===============================

function geo = makeSyntheticGeo(t0, dur_s)
% MAKESYNTHETICGEO  Fabricate a plausible overhead-pass geometry (selftest only).
%
% Inputs:
%   t0    - pass start (datetime, UTC)
%   dur_s - pass duration (s)
%
% Outputs:
%   geo - struct matching computePassGeometry's output shape:
%           .time (datetime UTC, 1 Hz), .az_deg, .el_deg (deg),
%           .range_m (m), .doppler_Hz (Hz at 435 MHz)
%
% Elevation follows a half-sine from -5 deg up to ~60 deg and back;
% Doppler follows the classic S-curve (+10 kHz -> -10 kHz through TCA).
t   = (0:1:dur_s).';
tau = t / dur_s;                                   % 0..1 through the pass
geo = struct();
geo.time       = t0 + seconds(t);
geo.el_deg     = -5 + 65 * sin(pi * tau);
geo.az_deg     = mod(180 + 180 * tau, 360);
geo.range_m    = 550e3 ./ sind(max(geo.el_deg, 7));   % crude slant range
geo.doppler_Hz = 10e3 * cos(pi * tau);
end
