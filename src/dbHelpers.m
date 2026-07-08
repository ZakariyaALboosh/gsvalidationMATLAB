function [db, undb] = dbHelpers()
% DBHELPERS  Single definition point for the dB conversion helpers.
%
%   The whole pipeline expresses power and SNR in dB. These two handles
%   are defined here and ONLY here; every module that needs a linear<->dB
%   conversion calls [db, undb] = dbHelpers() rather than redefining them.
%
% Inputs:
%   (none)
%
% Outputs:
%   db   - function handle, db(x)   = 10*log10(x)   (linear -> dB)
%   undb - function handle, undb(x) = 10.^(x/10)    (dB -> linear)

db   = @(x) 10*log10(x);
undb = @(x) 10.^(x/10);
end
