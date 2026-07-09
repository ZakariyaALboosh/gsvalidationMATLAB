function wf = parseWaterfallHDF5(path)
% PARSEWATERFALLHDF5  Read a SatNOGS network waterfall artifact (.h5/.hdf5).
%
%   Discovers the layout with h5info rather than assuming names, then maps
%   to the unified waterfall struct. Verified against a real artifact
%   (artifact_version 2), whose actual layout is:
%
%     /  attr 'metadata'  JSON: observation_id, frequency (STRING, Hz),
%                         tle (3 lines, '0 NAME' 3LE style), location
%     /waterfall          attrs: start_time, offset_in_stds, scale_in_stds
%       data          uint8, time x freq, quantised power in 'div' units
%       relative_time float64, seconds since start_time
%       frequency     float64, bin offsets from centre - the unit attribute
%                     claims kHz but the VALUES are Hz (checked: span
%                     +/-24000 with 46.875 spacing = 48 kHz / 1024 bins)
%       offset        float32, PER-CHANNEL dequantisation offset (dB)
%       scale         float32, PER-CHANNEL dequantisation scale (dB/div)
%
%   DEQUANTISATION IS MANDATORY: P_dB = offset(ch) + data*scale(ch).
%   The offset/scale vary per frequency channel (offset spans ~30 dB across
%   the passband), so the S - N differencing in extractMeasuredSNR would
%   NOT cancel them - only a global constant cancels. After dequantisation
%   the power is uncalibrated dB (arbitrary reference), which is exactly
%   what the SNR-only pipeline needs; no absolute calibration is attempted.
%
% Inputs:
%   path - path to the .h5/.hdf5 artifact
%
% Outputs:
%   wf - unified waterfall struct:
%          .t_s        Nx1 seconds since first row
%          .f_Hz       1xM bin offsets from centre (Hz)
%          .center_Hz  scalar centre frequency (Hz)
%          .P          NxM power (dB, uncalibrated)
%          .start_time datetime UTC of first row (NaT if absent)
%          .meta       struct: tle (3x1 cellstr), frequency_Hz, lat, lon,
%                      alt, observation_id

if exist(path, 'file') ~= 2
    error('parseWaterfallHDF5:missingFile', 'File not found: %s', path);
end

info = h5info(path);

% ---- root metadata JSON ------------------------------------------------
meta = struct();
if any(strcmp({info.Attributes.Name}, 'metadata'))
    md = jsondecode(char(h5readatt(path, '/', 'metadata')));
    if isfield(md, 'frequency')          % stored as a JSON *string*
        meta.frequency_Hz = str2double(string(md.frequency));
    end
    if isfield(md, 'tle')
        lines = strsplit(strtrim(char(md.tle)), {'\r\n', '\n'});
        meta.tle = lines(:);             % 3x1 cellstr, '0 NAME' kept verbatim
    end
    if isfield(md, 'location')
        meta.lat = md.location.latitude;
        meta.lon = md.location.longitude;
        meta.alt = md.location.altitude;
    end
    if isfield(md, 'observation_id')
        meta.observation_id = md.observation_id;
    end
else
    warning('parseWaterfallHDF5:noMetadata', ...
        'No root ''metadata'' attribute: TLE/location must come from config.');
end

% ---- locate the waterfall group ----------------------------------------
gidx = find(contains(lower({info.Groups.Name}), 'waterfall'), 1);
if isempty(gidx)
    error('parseWaterfallHDF5:noWaterfallGroup', ...
        'No group containing ''waterfall'' in %s. Groups found: %s', ...
        path, strjoin({info.Groups.Name}, ', '));
end
grp = info.Groups(gidx);
ds_names = {grp.Datasets.Name};

% ---- locate datasets: by expected name first, by shape as fallback ------
% spectrum: named 'data', else the only 2-D dataset in the group
spec_name = pickDataset(grp, ds_names, 'data', 2);
% time axis: named 'relative_time', else 1-D matching the long dim
t_name    = pickDataset(grp, ds_names, 'relative_time', 1);
% frequency axis: named 'frequency', else remaining 1-D
f_name    = pickDataset(grp, ds_names, 'frequency', 1);

P    = double(h5read(path, [grp.Name '/' spec_name]));
t_s  = double(h5read(path, [grp.Name '/' t_name]));
fax  = double(h5read(path, [grp.Name '/' f_name]));
t_s  = t_s(:);
fax  = fax(:).';

% ---- orient the spectrum to time x freq ---------------------------------
% The artifact stores (time x freq) in C order; MATLAB's h5read returns it
% transposed. Decide by matching dimensions against the axis lengths.
if size(P, 1) == numel(fax) && size(P, 2) == numel(t_s)
    P = P.';
elseif ~(size(P, 1) == numel(t_s) && size(P, 2) == numel(fax))
    error('parseWaterfallHDF5:shapeMismatch', ...
        'Spectrum is %dx%d but time axis has %d and freq axis %d samples.', ...
        size(P, 1), size(P, 2), numel(t_s), numel(fax));
end

% ---- dequantise (uint8 'div' units -> uncalibrated dB) ------------------
has_off = any(strcmp(ds_names, 'offset'));
has_scl = any(strcmp(ds_names, 'scale'));
if has_off && has_scl
    off = double(h5read(path, [grp.Name '/offset'])); off = off(:).';
    scl = double(h5read(path, [grp.Name '/scale']));  scl = scl(:).';
    if numel(off) ~= size(P, 2)
        error('parseWaterfallHDF5:badQuantAxes', ...
            'offset/scale have %d entries but the spectrum has %d channels.', ...
            numel(off), size(P, 2));
    end
    P = off + P .* scl;              % per-channel: see header, mandatory
elseif max(P(:)) <= 255 && min(P(:)) >= 0
    error('parseWaterfallHDF5:quantisedWithoutScale', ...
        ['Spectrum looks uint8-quantised but the group has no offset/', ...
         'scale datasets - cannot reconstruct dB. Datasets found: %s'], ...
        strjoin(ds_names, ', '));
end
% (float-valued artifacts without offset/scale pass through unchanged)

% ---- frequency axis units -----------------------------------------------
% The unit attribute claims kHz but real files carry Hz. Heuristic on the
% span: a SatNOGS waterfall is tens of kHz wide, so a max |value| under
% 1000 can only be kHz; anything larger is already Hz.
if max(abs(fax)) < 1000
    fax = fax * 1e3;
end

% ---- start_time attribute ------------------------------------------------
start_time = NaT('TimeZone', 'UTC');
attr_names = {grp.Attributes.Name};
if any(strcmp(attr_names, 'start_time'))
    raw = strtrim(char(h5readatt(path, grp.Name, 'start_time')));
    start_time = parseISOTime(raw);
end

% ---- assemble --------------------------------------------------------------
wf = struct();
wf.t_s        = t_s - t_s(1);
wf.f_Hz       = fax;
wf.center_Hz  = NaN;
if isfield(meta, 'frequency_Hz'), wf.center_Hz = meta.frequency_Hz; end
wf.P          = P;
wf.start_time = start_time;
wf.meta       = meta;
end

function name = pickDataset(grp, ds_names, wanted, ndims_wanted)
% PICKDATASET  Dataset by expected name, else by dimensionality (defensive).
% Inputs:  grp (h5info group), ds_names (cellstr), wanted (char),
%          ndims_wanted (1 or 2 - axis vs spectrum)
% Outputs: name (char) - errors if nothing matches
if any(strcmp(ds_names, wanted))
    name = wanted;
    return;
end
sizes = {grp.Datasets.Dataspace};
is_nd = cellfun(@(s) numel(s.Size(s.Size > 1)) == ndims_wanted, sizes);
idx = find(is_nd, 1);
if isempty(idx)
    error('parseWaterfallHDF5:datasetNotFound', ...
        'No dataset named ''%s'' and no %d-D fallback. Found: %s', ...
        wanted, ndims_wanted, strjoin(ds_names, ', '));
end
name = ds_names{idx};
warning('parseWaterfallHDF5:nameFallback', ...
    'Dataset ''%s'' not found; using ''%s'' by shape.', wanted, name);
end

function t = parseISOTime(raw)
% PARSEISOTIME  ISO-8601 UTC string (with or without fraction) -> datetime UTC.
% Inputs:  raw - e.g. '2026-06-11T18:24:14.782770Z'
% Outputs: t   - datetime with TimeZone UTC; NaT if unparseable
fmts = {'yyyy-MM-dd''T''HH:mm:ss.SSSSSS''Z''', ...
        'yyyy-MM-dd''T''HH:mm:ss''Z''', ...
        'yyyy-MM-dd''T''HH:mm:ss.SSS''Z'''};
for k = 1:numel(fmts)
    try                              % datetime THROWS on format mismatch
        t = datetime(raw, 'InputFormat', fmts{k}, 'TimeZone', 'UTC');
        if ~isnat(t), return; end
    catch
    end
end
warning('parseWaterfallHDF5:badStartTime', ...
    'Could not parse start_time ''%s''; returning NaT.', raw);
t = NaT('TimeZone', 'UTC');
end
