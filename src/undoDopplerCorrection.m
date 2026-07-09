function wf_out = undoDopplerCorrection(wf, geo, opts)
% UNDODOPPLERCORRECTION  Reconstruct the raw sky-frequency waterfall.
%
%   SatNOGS waterfalls are recorded AFTER the client's Doppler correction:
%   the displayed frequency of every row is shifted by -doppler(t), so the
%   satellite appears vertical and fixed carriers draw inverted S-curves.
%   This utility undoes that shift row by row, recovering the waterfall as
%   it would look at the antenna: the satellite becomes the classic
%   Doppler S-curve and fixed carriers become vertical lines.
%
%   Mapping (see extractMeasuredSNR header for the empirical basis):
%     corrected frequency  nu = f_true - doppler(t)
%     =>  P_raw(t, f) = P_corr(t, f - doppler(t))     [row-wise resample]
%
%   The resample is linear interpolation on the fixed bin grid. Bins that
%   shift in from outside the recorded passband are unknowable; they are
%   filled with the row's median (a flat noise estimate) so downstream
%   statistics stay usable - do not treat the outer +/-|doppler| edge
%   strips as real data.
%
% Inputs:
%   wf   - Doppler-CORRECTED waterfall struct from parseWaterfall:
%            .t_s (Nx1 s), .f_Hz (1xM Hz), .P (NxM dB), .start_time
%   geo  - geometry struct from computePassGeometry (drives doppler(t)):
%            .time (datetime UTC), .doppler_Hz (Hz)
%   opts - optional struct, any subset of:
%            .time_offset_s manual waterfall->geometry alignment nudge (s),
%                           same meaning as in extractMeasuredSNR [default 0]
%            .out_h5        path to also write the result as a simple HDF5
%                           file (float data + axes + metadata; readable by
%                           parseWaterfallHDF5)            [default: no file]
%
% Outputs:
%   wf_out - same struct layout as wf, with .P de-corrected. Rows outside
%            the geometry time span (doppler unknown) are copied unchanged
%            and counted in a warning.

if nargin < 3, opts = struct(); end
if ~isfield(opts, 'time_offset_s'), opts.time_offset_s = 0;  end
if ~isfield(opts, 'out_h5'),        opts.out_h5        = ''; end

% ---- doppler on the waterfall row times (same alignment as extraction) --
geo_t_s = seconds(geo.time(:) - geo.time(1));
if isfield(wf, 'start_time') && isdatetime(wf.start_time) && ~isnat(wf.start_time)
    toff = seconds(wf.start_time - geo.time(1)) + opts.time_offset_s;
else
    toff = opts.time_offset_s;
end
dopp_Hz = interp1(geo_t_s, geo.doppler_Hz(:), wf.t_s(:) + toff, 'linear', NaN);

% ---- row-wise inverse shift ---------------------------------------------
wf_out = wf;
f_Hz   = wf.f_Hz(:).';
n_skipped = 0;
for i = 1:numel(wf.t_s)
    if ~isfinite(dopp_Hz(i))
        n_skipped = n_skipped + 1;      % no doppler -> cannot de-correct
        continue;
    end
    row = interp1(f_Hz, wf.P(i, :), f_Hz - dopp_Hz(i), 'linear', NaN);
    row(~isfinite(row)) = median(wf.P(i, :));   % edge fill: see header
    wf_out.P(i, :) = row;
end
if n_skipped > 0
    warning('undoDopplerCorrection:rowsOutsideGeometry', ...
        '%d of %d rows lie outside the geometry time span; left unchanged.', ...
        n_skipped, numel(wf.t_s));
end

% ---- optional HDF5 export ------------------------------------------------
if ~isempty(opts.out_h5)
    if exist(opts.out_h5, 'file') == 2, delete(opts.out_h5); end
    % Written in MATLAB (column-major) order; parseWaterfallHDF5 resolves
    % the orientation by matching axis lengths, so it reads back correctly.
    h5create(opts.out_h5, '/waterfall/data', size(wf_out.P), 'Datatype', 'single');
    h5write(opts.out_h5,  '/waterfall/data', single(wf_out.P));
    h5create(opts.out_h5, '/waterfall/relative_time', numel(wf_out.t_s));
    h5write(opts.out_h5,  '/waterfall/relative_time', wf_out.t_s(:));
    h5create(opts.out_h5, '/waterfall/frequency', numel(wf_out.f_Hz));
    h5write(opts.out_h5,  '/waterfall/frequency', wf_out.f_Hz(:));
    if isdatetime(wf_out.start_time) && ~isnat(wf_out.start_time)
        h5writeatt(opts.out_h5, '/waterfall', 'start_time', ...
            char(wf_out.start_time, 'yyyy-MM-dd''T''HH:mm:ss.SSSSSS''Z''')); %#ok<DATST>
    end
    if ~isempty(fieldnames(wf_out.meta))
        % re-emit in the SatNOGS artifact schema so parseWaterfallHDF5 can
        % read the exported file back (frequency/tle/location field names)
        md = struct();
        m  = wf_out.meta;
        if isfield(m, 'observation_id'), md.observation_id = m.observation_id; end
        if isfield(m, 'frequency_Hz'),   md.frequency = m.frequency_Hz;        end
        if isfield(m, 'tle'),            md.tle = strjoin(m.tle, newline);     end
        if isfield(m, 'lat')
            md.location = struct('latitude', m.lat, 'longitude', m.lon, ...
                                 'altitude', m.alt);
        end
        h5writeatt(opts.out_h5, '/', 'metadata', jsonencode(md));
    end
    fprintf('[undoDopplerCorrection] wrote de-corrected waterfall to %s\n', ...
        opts.out_h5);
end
end
