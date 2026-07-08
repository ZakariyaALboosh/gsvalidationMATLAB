function [snr_pred_dB, budget, excluded] = predictSNR(geo, station)
% PREDICTSNR  Link-budget SNR prediction over a pass.
%
%   Per time step:
%     FSPL_dB    = 20*log10(4*pi*range*f/c)
%     L_atm_dB   = gaspl (standard atmosphere, cosecant path) if available
%                  and f >= 1 GHz (gaspl's validity floor); else a fixed
%                  0.2 dB with a warning
%     Prx_dBW    = EIRP - FSPL - L_atm - line_loss - pol_loss + G(el)
%     N_dBW      = 10*log10(k * T_sys * B)
%     SNR_dB     = Prx - N
%
% Inputs:
%   geo     - struct from computePassGeometry:
%               .el_deg (deg), .range_m (m), .time (datetime UTC)
%   station - struct from loadStationConfig:
%               .freq_Hz (Hz), .sat_eirp_dBW (dBW), .ant_gain_dBi (dBi),
%               .ant_pattern (@(el_deg)->dBi), .sys_noise_temp_K (K),
%               .rx_bw_Hz (Hz), .line_loss_dB (dB), .pol_loss_dB (dB)
%
% Outputs:
%   snr_pred_dB - predicted SNR aligned to geo.time (dB, Nx1)
%   budget      - struct of the constant budget terms (for the thesis
%                 table): freq_Hz, sat_eirp_dBW, ant_gain_dBi,
%                 line_loss_dB, pol_loss_dB, sys_noise_temp_K, rx_bw_Hz,
%                 noise_dBW, atm_model, fspl_min_dB, fspl_max_dB
%   excluded    - logical Nx1, true where el < 5 deg (do not use these
%                 samples in the comparison)

[db, ~] = dbHelpers();
c = 299792458;                 % m/s
k = 1.380649e-23;              % J/K

el = geo.el_deg(:);
f  = station.freq_Hz;

% ---- path terms --------------------------------------------------------
fspl_dB = 20 * log10(4 * pi * geo.range_m(:) * f / c);

if exist('gaspl', 'file') == 2 && f >= 1e9
    % standard atmosphere, path length ~ atmosphere thickness / sin(el)
    path_m   = 8000 ./ sind(max(el, 5));
    L_atm_dB = gaspl(path_m, f, 288.15, 101325, 7.5);
    atm_model = 'gaspl, std atmosphere, cosecant path';
else
    L_atm_dB = 0.2 * ones(size(el));
    atm_model = 'fixed 0.2 dB';
    warning('predictSNR:atmFallback', ...
        ['gaspl unavailable or f < 1 GHz (its validity floor): using a ', ...
         'fixed 0.2 dB atmospheric loss.']);
end

% ---- antenna gain vs elevation -----------------------------------------
if isfield(station, 'ant_pattern') && ~isempty(station.ant_pattern)
    G_dBi = station.ant_pattern(el);
else
    G_dBi = station.ant_gain_dBi * ones(size(el));
end

% ---- received power, noise, SNR ----------------------------------------
Prx_dBW = station.sat_eirp_dBW - fspl_dB - L_atm_dB ...
          - station.line_loss_dB - station.pol_loss_dB + G_dBi;
N_dBW   = db(k * station.sys_noise_temp_K * station.rx_bw_Hz);

snr_pred_dB = Prx_dBW - N_dBW;
excluded    = el < 5;

% ---- constant terms for the thesis table --------------------------------
budget = struct();
budget.freq_Hz          = f;
budget.sat_eirp_dBW     = station.sat_eirp_dBW;
budget.ant_gain_dBi     = station.ant_gain_dBi;
budget.line_loss_dB     = station.line_loss_dB;
budget.pol_loss_dB      = station.pol_loss_dB;
budget.sys_noise_temp_K = station.sys_noise_temp_K;
budget.rx_bw_Hz         = station.rx_bw_Hz;
budget.noise_dBW        = N_dBW;
budget.atm_model        = atm_model;
budget.fspl_min_dB      = min(fspl_dB);
budget.fspl_max_dB      = max(fspl_dB);
end
