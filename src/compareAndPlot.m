function results = compareAndPlot(geo, snr_pred_dB, wf, snr_meas_dB, info)
% COMPAREANDPLOT  Offset-fit comparison of predicted vs measured SNR + figures.
%
%   Interpolates the prediction onto the measurement timestamps, restricts
%   to elevation > min_el_deg, and computes:
%     - mean offset      = mean(meas - pred)   (dB)  [the fitted station
%                          constant: absorbs unknown EIRP / G/T / scaling]
%     - RMS residual     after removing that offset (dB)
%     - Pearson correlation of the two curves
%   Saves three figures (PNG + .fig) and a results .mat to the output dir,
%   and prints a one-paragraph console interpretation.
%
% Inputs:
%   geo         - geometry struct: .time (datetime UTC), .el_deg (deg)
%   snr_pred_dB - predicted SNR on geo.time (dB, same length as geo.time)
%   wf          - waterfall struct (.t_s s, .f_Hz Hz, .P dB, .start_time)
%   snr_meas_dB - measured SNR on wf.t_s (dB, from extractMeasuredSNR)
%   info        - optional struct, any subset of:
%                   .out_dir       output directory   [default '../output'
%                                  relative to this file's folder]
%                   .label         tag used in filenames/titles ['run']
%                   .time_offset_s effective waterfall->geometry offset (s)
%                                  as returned by extractMeasuredSNR [auto:
%                                  recomputed from wf.start_time, else 0]
%                   .track         Nx1 signal-track bin indices (for the
%                                  waterfall overlay figure)      [none]
%                   .min_el_deg    elevation mask (deg)           [5]
%
% Outputs:
%   results - struct:
%               .offset_dB   mean(meas - pred) over the mask (dB)
%               .rms_dB      RMS of (meas - pred - offset) (dB)
%               .corr        Pearson correlation coefficient
%               .n_used      number of samples in the mask
%               .t_meas_s    measurement times on the geometry axis (s)
%               .snr_meas_dB, .snr_pred_i_dB, .el_i_deg  per-sample vectors
%               .mask        logical vector of samples used
%               .residual_dB residuals after offset removal (masked samples
%                            only, NaN elsewhere)

if nargin < 5, info = struct(); end
if ~isfield(info, 'out_dir') || isempty(info.out_dir)
    info.out_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'output');
end
if ~isfield(info, 'label'),      info.label      = 'run'; end
if ~isfield(info, 'min_el_deg'), info.min_el_deg = 5;     end
if ~isfield(info, 'time_offset_s') || isempty(info.time_offset_s)
    if isfield(wf, 'start_time') && isdatetime(wf.start_time) && ~isnat(wf.start_time)
        info.time_offset_s = seconds(wf.start_time - geo.time(1));
    else
        info.time_offset_s = 0;
    end
end
if ~exist(info.out_dir, 'dir'), mkdir(info.out_dir); end

% ---- put measurement rows on the geometry time axis -------------------
geo_t_s = seconds(geo.time(:) - geo.time(1));
t_meas  = wf.t_s(:) + info.time_offset_s;
snr_meas_dB = snr_meas_dB(:);

snr_pred_i = interp1(geo_t_s, snr_pred_dB(:), t_meas, 'linear', NaN);
el_i       = interp1(geo_t_s, geo.el_deg(:),  t_meas, 'linear', NaN);

mask = isfinite(snr_meas_dB) & isfinite(snr_pred_i) & (el_i > info.min_el_deg);
if nnz(mask) < 2
    error('compareAndPlot:noOverlap', ...
        ['Fewer than 2 valid samples above %g deg elevation. Check the ', ...
         'time alignment (time_offset_s = %.1f s) and the pass window.'], ...
        info.min_el_deg, info.time_offset_s);
end

% ---- offset-fit statistics --------------------------------------------
d          = snr_meas_dB - snr_pred_i;                 % meas - pred (dB)
offset_dB  = mean(d(mask));
residual   = nan(size(d));
residual(mask) = d(mask) - offset_dB;
rms_dB     = sqrt(mean(residual(mask).^2));
cc         = corrcoef(snr_meas_dB(mask), snr_pred_i(mask));
corr_val   = cc(1, 2);

results = struct();            % field-by-field: struct(...) with vector
results.offset_dB     = offset_dB;   % values would build a struct ARRAY
results.rms_dB        = rms_dB;
results.corr          = corr_val;
results.n_used        = nnz(mask);
results.t_meas_s      = t_meas;
results.snr_meas_dB   = snr_meas_dB;
results.snr_pred_i_dB = snr_pred_i;
results.el_i_deg      = el_i;
results.mask          = mask;
results.residual_dB   = residual;

% ---- figure 1: SNR curves vs time, elevation on right axis ------------
fig1 = figure('Visible', 'off', 'Name', 'SNR vs time');
yyaxis left
plot(geo_t_s, snr_pred_dB(:), '-', 'LineWidth', 1.2); hold on
plot(t_meas(mask), snr_meas_dB(mask), '.', 'MarkerSize', 8);
ylabel('SNR (dB)');
yyaxis right
plot(geo_t_s, geo.el_deg(:), '--');
ylabel('Elevation (deg)');
xlabel('Time since pass start (s)');
title(sprintf('%s: predicted vs measured SNR', info.label), 'Interpreter', 'none');
legend({'predicted', 'measured', 'elevation'}, 'Location', 'south');
grid on
saveFig(fig1, info.out_dir, [info.label '_snr_vs_time']);

% ---- figure 2: residuals vs elevation ----------------------------------
fig2 = figure('Visible', 'off', 'Name', 'Residuals vs elevation');
plot(el_i(mask), residual(mask), '.', 'MarkerSize', 8); hold on
yline(0, '-');
xlabel('Elevation (deg)');
ylabel(sprintf('Residual after %.2f dB offset removal (dB)', offset_dB));
title(sprintf('%s: residuals vs elevation (RMS %.2f dB)', info.label, rms_dB), ...
    'Interpreter', 'none');
grid on
saveFig(fig2, info.out_dir, [info.label '_residual_vs_el']);

% ---- figure 3: waterfall with Doppler track overlay --------------------
fig3 = figure('Visible', 'off', 'Name', 'Waterfall + track');
imagesc(wf.f_Hz / 1e3, wf.t_s, wf.P); hold on
axis xy
if isfield(info, 'track') && ~isempty(info.track)
    tr = info.track(:);
    ok = isfinite(tr);
    plot(wf.f_Hz(tr(ok)) / 1e3, wf.t_s(ok), 'r-', 'LineWidth', 1.0);
end
xlabel('Frequency offset from centre (kHz)');
ylabel('Time since first row (s)');
title(sprintf('%s: waterfall with predicted Doppler track', info.label), ...
    'Interpreter', 'none');
cb = colorbar; cb.Label.String = 'Power (dB, uncalibrated)';
saveFig(fig3, info.out_dir, [info.label '_waterfall_track']);

% ---- results .mat + console summary ------------------------------------
save(fullfile(info.out_dir, [info.label '_results.mat']), 'results', 'info');

if abs(offset_dB) < 1 && rms_dB < 2
    interp_txt = ['offset ~0 with a small, structureless residual: the ', ...
        'link budget is validated in an absolute sense for this station.'];
else
    interp_txt = sprintf(['the nonzero offset is reported as the fitted ', ...
        'station constant (%.2f dB): it absorbs the unknown EIRP / G/T / ', ...
        'waterfall scaling; the SHAPE agreement (RMS, correlation) is the ', ...
        'validation result.'], offset_dB);
end
fprintf(['\n[compareAndPlot] %s: %d samples above %g deg elevation. ', ...
    'Mean offset (meas-pred) = %.2f dB; RMS residual after offset ', ...
    'removal = %.2f dB; Pearson correlation = %.3f. Interpretation: %s\n'], ...
    info.label, nnz(mask), info.min_el_deg, offset_dB, rms_dB, corr_val, ...
    interp_txt);
end

function saveFig(fig, out_dir, name)
% SAVEFIG  Save a figure as PNG and .fig into out_dir (helper, no units).
print(fig, fullfile(out_dir, [name '.png']), '-dpng', '-r150');
savefig(fig, fullfile(out_dir, [name '.fig']));
close(fig);
end
