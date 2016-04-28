function [ IDstr, x, xnorm, y, p, cp, ref ] = Load_NACA_Data( path, data_filename )
%%%
% Loads data from a CFD run saved in ParaView, and calculates Cp, PR, etc.
%
% Usage:
%   Load_CFD_Data('\directory\of\files\','specific_file.csv')
%    *ASSUMES path is backslash-terminated
%%%

    % Assert that inputs are nonempty strings.
    validateattributes(         path,{'char'},{'nonempty'});
    validateattributes(data_filename,{'char'},{'nonempty'});
    
    %%%
    % Determine full path to data file and reference file (upstream pressure probe).
    %%%
    
    data_path = [path, data_filename];
    
    under_indices = strfind(data_filename,'_');
    common_prefix = data_filename(1:under_indices(end));
    ref_path = [path, common_prefix, 'probe.csv'];
    
    IDstr = data_filename(under_indices(1)+1:under_indices(2)-1);
    
    %%%
    % Load raw data from ParaView files.
    %%%
    
    [x, y, p] = Load_ParaView_Data(data_path, {'coordsX', 'coordsY', 'pressure'});
    
    [ref.x, ref.y, ref.p, ref.u] = ...
                Load_ParaView_Data(ref_path,  {'coordsX', 'coordsY', 'pressure', ...
                                               'velocity:0'});
    
	%%%
    % Compute extra fields relevant to flow analysis.
    %%%
    
    % Normalized x-coordinate
    xnorm = x - min(x);
    xnorm = xnorm / max(xnorm);
    
    % Physical constants.
	gamma = 1.4;
    rho = 1.0;
    R = 288.294;
    
    % Flow pressure coefficient.
    cp = (p - ref.p) ./ (0.5 * rho * (ref.u)^2);
    
%     % Flow pressure recovery.
%     uMag = sqrt(u.^2 + v.^2 + w.^2);
%     numerator = p .* (1 + 0.5 * uMag.^2 ./ (R * T));
%     denomenator = ref.p * (1 + 0.5 * norm(ref.u)^2 / (R * ref.T));
%     pr = numerator / denomenator;

end