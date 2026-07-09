function wf = parseWaterfall(path)
% PARSEWATERFALL  Dispatch a waterfall file to the parser for its format.
%
%   .h5 / .hdf5  -> parseWaterfallHDF5 (SatNOGS network artifact)
%   .dat         -> parseWaterfallDAT  (raw client file; build step 4)
%
% Inputs:
%   path - path to the waterfall file
%
% Outputs:
%   wf - unified waterfall struct (see parseWaterfallHDF5 header):
%          .t_s (Nx1 s), .f_Hz (1xM Hz), .center_Hz (Hz), .P (NxM dB),
%          .start_time (datetime UTC or NaT), .meta (struct; empty for .dat)

[~, ~, ext] = fileparts(path);
switch lower(ext)
    case {'.h5', '.hdf5'}
        wf = parseWaterfallHDF5(path);
    case '.dat'
        if exist('parseWaterfallDAT', 'file') ~= 2
            error('parseWaterfall:notBuilt', ...
                'parseWaterfallDAT arrives with build step 4.');
        end
        wf = parseWaterfallDAT(path);
    otherwise
        error('parseWaterfall:badExtension', ...
            'Unsupported waterfall extension ''%s'' (want .h5/.hdf5/.dat).', ext);
end
end
