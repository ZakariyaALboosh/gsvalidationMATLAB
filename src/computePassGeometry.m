function geo = computePassGeometry(tle, lat_deg, lon_deg, alt_m, t_start, t_stop, step_s, freq_Hz)
% COMPUTEPASSGEOMETRY  Az/el/range/Doppler histories for one station and pass.
%
%   Builds a satelliteScenario (Aerospace Toolbox / Satellite Communications
%   Toolbox), propagates the TLE, and returns the ground-station-relative
%   geometry over [t_start, t_stop]. Doppler is computed from the numerical
%   range rate.
%
% Inputs:
%   tle      - TLE as a file path (char/string), OR the TLE text itself
%              (3-line char matrix, 3-element string array, or one char/
%              string containing newlines). Text input is written to a
%              temp file because satellite() reads TLEs from file.
%   lat_deg  - station geodetic latitude (deg, +N)
%   lon_deg  - station longitude (deg, +E)
%   alt_m    - station altitude above WGS84 ellipsoid (m)
%   t_start  - observation start (datetime, UTC)
%   t_stop   - observation stop (datetime, UTC)
%   step_s   - sample step (s)                                  [default 1]
%   freq_Hz  - downlink frequency for the Doppler conversion (Hz)
%
% Outputs:
%   geo - struct:
%           .time       datetime vector, UTC (Nx1)
%           .az_deg     azimuth (deg, Nx1)
%           .el_deg     elevation (deg, Nx1)
%           .range_m    slant range (m, Nx1)
%           .doppler_Hz Doppler shift at freq_Hz (Hz, Nx1; + = approaching)

if nargin < 7 || isempty(step_s), step_s = 1; end

if exist('satelliteScenario', 'file') ~= 2
    error('computePassGeometry:missingToolbox', ...
        ['satelliteScenario not found. This module needs the Aerospace ', ...
         'Toolbox or Satellite Communications Toolbox.']);
end

% ---- normalise the TLE to a file (satellite() reads TLEs from file) ---
tle_file = resolveTLEFile(tle);

% ---- scenario times: force UTC wall-clock, then strip the zone --------
% (satelliteScenario treats unzoned datetimes as UTC; passing the zone
% explicitly is version-sensitive, so convert first.)
t_start.TimeZone = 'UTC';  t_start.TimeZone = '';
t_stop.TimeZone  = 'UTC';  t_stop.TimeZone  = '';

sc  = satelliteScenario(t_start, t_stop, step_s);
sat = satellite(sc, tle_file);
gs  = groundStation(sc, lat_deg, lon_deg, 'Altitude', alt_m);

% aer with no time argument returns full histories over the scenario
% timeline (AutoSimulate default); timeOut is the matching datetime vector.
[az, el, rng_m, t_out] = aer(gs, sat);

geo = struct();
geo.time    = t_out(:);
geo.time.TimeZone = 'UTC';
geo.az_deg  = az(:);
geo.el_deg  = el(:);
geo.range_m = rng_m(:);

if max(geo.el_deg) < 0
    error('computePassGeometry:noPass', ...
        ['Satellite never rises above the horizon in %s .. %s (max el ', ...
         '%.1f deg). Check that the TLE matches the observation and that ', ...
         'the time window actually contains the pass.'], ...
        string(t_start), string(t_stop), max(geo.el_deg));
end

% ---- Doppler from numerical range rate --------------------------------
c  = 299792458;                                            % m/s
rr = gradient(geo.range_m) ./ gradient(seconds(geo.time - geo.time(1)));
geo.doppler_Hz = -rr / c * freq_Hz;
end

function tle_file = resolveTLEFile(tle)
% RESOLVETLEFILE  Return a readable TLE file path for any accepted TLE input.
%
% Inputs:  tle - file path, or TLE text (see computePassGeometry header)
% Outputs: tle_file - char path to an existing TLE file
if (ischar(tle) && size(tle, 1) == 1 || (isstring(tle) && isscalar(tle))) ...
        && exist(tle, 'file') == 2
    tle_file = char(tle);
    return;
end
% Not a file: treat as TLE text and normalise to a 3-line cellstr
if iscellstr(tle) && numel(tle) == 3                       %#ok<ISCLSTR>
    lines = tle(:).';
elseif isstring(tle) && numel(tle) == 3
    lines = cellstr(tle(:).');
elseif ischar(tle) && size(tle, 1) == 3
    lines = cellstr(tle);
elseif (ischar(tle) || (isstring(tle) && isscalar(tle)))
    lines = strsplit(strtrim(char(tle)), {'\r\n', '\n'});
else
    error('computePassGeometry:badTLE', ...
        'TLE must be a file path or 3 lines of TLE text.');
end
if numel(lines) ~= 3
    error('computePassGeometry:badTLE', ...
        'TLE text has %d lines; expected 3 (name + line1 + line2).', ...
        numel(lines));
end
% SatNOGS metadata uses the 3LE name form '0 SWISSCUBE'; strip the '0 '
% marker so the name line is plain (safest for MATLAB's TLE reader)
if startsWith(lines{1}, '0 ')
    lines{1} = strtrim(lines{1}(3:end));
end
tle_file = [tempname '.tle'];
fid = fopen(tle_file, 'w');
fprintf(fid, '%s\n', lines{:});
fclose(fid);
end
