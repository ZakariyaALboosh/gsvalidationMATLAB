function station = loadStationConfig(config_path)
% LOADSTATIONCONFIG  Load and validate the station configuration script.
%
%   Runs the config script (which must fill a struct named `station`),
%   checks required fields, and applies defaults.
%
%   PRECEDENCE (documented also in the config file): when a waterfall's
%   wf.meta is populated (SatNOGS HDF5 artifact), the artifact's location,
%   frequency and TLE OVERRIDE this config for the geometry. The config
%   always supplies the RF terms: ant_gain_dBi / ant_pattern,
%   sys_noise_temp_K, rx_bw_Hz, line_loss_dB, pol_loss_dB, sat_eirp_dBW.
%
% Inputs:
%   config_path - path to the config script, e.g. ../config/station_uhf.m
%
% Outputs:
%   station - struct with fields:
%               lat_deg (deg), lon_deg (deg), alt_m (m), freq_Hz (Hz),
%               ant_gain_dBi (dBi), ant_pattern (@(el_deg)->dBi or []),
%               sys_noise_temp_K (K), rx_bw_Hz (Hz), line_loss_dB (dB),
%               pol_loss_dB (dB), sat_eirp_dBW (dBW)

if exist(config_path, 'file') ~= 2
    error('loadStationConfig:missingFile', ...
        'Station config not found: %s', config_path);
end

run(config_path);   % config script defines `station` in this workspace

if ~exist('station', 'var') || ~isstruct(station)
    error('loadStationConfig:badConfig', ...
        '%s must define a struct named ''station''.', config_path);
end

% ---- defaults -----------------------------------------------------------
if ~isfield(station, 'freq_Hz') || isempty(station.freq_Hz)
    station.freq_Hz = 435e6;
end
if ~isfield(station, 'ant_pattern')
    station.ant_pattern = [];
end

% ---- required fields ----------------------------------------------------
required = {'lat_deg', 'lon_deg', 'alt_m', 'ant_gain_dBi', ...
    'sys_noise_temp_K', 'rx_bw_Hz', 'line_loss_dB', 'pol_loss_dB', ...
    'sat_eirp_dBW'};
missing = required(~isfield(station, required));
if ~isempty(missing)
    error('loadStationConfig:missingFields', ...
        'Station config %s is missing field(s): %s', ...
        config_path, strjoin(missing, ', '));
end
end
