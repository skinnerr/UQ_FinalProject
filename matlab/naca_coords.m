function [yt, xu, yu, xl, yl] = naca_coords(x, m, p, t)
%%%
% Generates upper (xu, yu) and lower (xl, yl) surfaces of a 4-digit NACA airfoil.
%
%   x = coordinates along the chord (x-axis), s.t. max(x) = chord length.
%   m = max camber
%   p = location of max camber over chord
%   t = maximum thickness over chord
%%%
    
    c = max(x);
    ipc = find(x > p*c, 1); % First index where x > p*c
    
    % Half-thickness for an un-cambered airfoil
    yt = 5*t*c * ( 0.2969*(x/c).^0.5 + ...
                  -0.1260*(x/c)      + ...
                  -0.3516*(x/c).^2   + ...
                   0.2843*(x/c).^3   + ...
                  -0.1015*(x/c).^4     );

	% Camber line
    yc = nan(size(x));
    yc(1:ipc-1) = (m*(  x(1:ipc-1))/(  p)^2).*(  2*p-x(1:ipc-1)/c);
    yc(ipc:end) = (m*(c-x(ipc:end))/(1-p)^2).*(1-2*p+x(ipc:end)/c);
    
    % Derivative of the camber line y_c w.r.t. space
    dcamb = nan(size(x));
    dcamb(1:ipc-1) = (2*m/(  p)^2)*(p-x(1:ipc-1)/c);
    dcamb(ipc:end) = (2*m/(1-p)^2)*(p-x(ipc:end)/c);
    theta = atan(dcamb);
    
    % Coordinate pairs of upper and lower cambered surfaces (xu, yu) and (xl, yl)
    xu = x - yt .* sin(theta);
    xl = x + yt .* sin(theta);
    yu = yc + yt .* cos(theta);
    yl = yc - yt .* cos(theta);
    
end










