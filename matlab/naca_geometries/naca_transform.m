function [dx, dy] = naca_transform(x0, y0, m, p, t, c, a)
%%%
% Calculates the displacement vector to transform a point on a NACA 0012 airfoil of unit
% chord to an arbitrary 4-digit NACA airfoil at a given angle of attack.
%
%   x0 = x-coordinate of reference point
%   y0 = y-coordinate of reference point
%   m  = max camber
%   p  = location of max camber over chord
%   t  = maximum thickness over chord
%   c  = chord length
%   a  = angle of attack
%%%

    % Properties of the reference unit-chord NACA 0012 airfoil.
    c0 = 1;
    
    % Scale our actual x-coordinate to account for change in chord.
    x = x0 * c / c0;
    
    % Half-thickness describing symmetric airfoil.
    yt = 5*t*c * ( 0.2969*(x/c).^0.5 + ...
                  -0.1260*(x/c)      + ...
                  -0.3516*(x/c).^2   + ...
                   0.2843*(x/c).^3   + ...
                  -0.1015*(x/c).^4     );
    
    % Determine region of point relative to location of max camber.
    xgt = x > p*c;

	% Camber line and its derivative w.r.t. x.
    if xgt
        yc = (m*(c-x)/(1-p)^2)*(1-2*p+x/c);
        dc = (2*m/(1-p)^2)*(p-x/c);
    else
        yc = (m*(  x)/(  p)^2)*(  2*p-x/c);
        dc = (2*m/(  p)^2)*(p-x/c);
    end
    th = atan(dc);
    
    % Coordinate pairs on new surface, using the sign of y0 to determine initial surface.
    if sign(y0) >= 0
        x = x - yt * sin(th);
        y = yc + yt * cos(th);
    else
        x = x + yt * sin(th);
        y = yc - yt * cos(th);
    end
    
    % Change sign of a so it specifies a clockwise rotation of the airfoil, convert to
    % radians.
    a = -a * pi / 180;
    
    % Rotate new surface coordinate pairs for angle of attack.
    rotmat = [cos(a), -sin(a); sin(a), cos(a)];
    tmp = rotmat * [x;y];
    x = tmp(1);
    y = tmp(2);
    
    
    dx = x - x0;
    dy = y - y0;
    
end










