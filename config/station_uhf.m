% STATION_UHF  Station configuration (edit the placeholder values).
%
%   Plain script: fills a struct named `station`, read by loadStationConfig.
%
%   PRECEDENCE: when the waterfall is a SatNOGS HDF5 artifact, its embedded
%   metadata (location, frequency, TLE) OVERRIDES the values here for the
%   GEOMETRY. This file always supplies the RF terms (gains, losses, noise
%   temperature, EIRP), which the artifact does not carry.

station = struct();

% ---- location (used only when the waterfall has no metadata, e.g. .dat) --
station.lat_deg = 52.0;        % EDIT ME: geodetic latitude (deg, +N)
station.lon_deg = 4.4;         % EDIT ME: longitude (deg, +E)
station.alt_m   = 10;          % EDIT ME: altitude (m)

% ---- link ---------------------------------------------------------------
station.freq_Hz = 435e6;       % downlink centre frequency (Hz)

% ---- RF terms (always taken from this file) ------------------------------
station.ant_gain_dBi     = 12;   % EDIT ME: antenna boresight gain (dBi)
station.ant_pattern      = [];   % optional @(el_deg)->dBi; [] = constant gain
station.sys_noise_temp_K = 500;  % EDIT ME: system noise temperature (K)
station.rx_bw_Hz         = 5000; % noise bandwidth for the SNR (Hz); keep
                                 % equal to sig_bw_Hz in extractMeasuredSNR
station.line_loss_dB     = 1;    % EDIT ME: feedline/connector loss (dB)
station.pol_loss_dB      = 3;    % EDIT ME: polarisation mismatch loss (dB)
station.sat_eirp_dBW     = -3;   % EDIT ME: satellite EIRP (dBW)
