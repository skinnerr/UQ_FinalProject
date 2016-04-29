function [ camber_line ] = NACA_Camber_Line( x, m, p, c0, c, a )
%%%
% Calculates the camber line y-coordinates of an airfoil with the given parameters.
%
%   x  = coordinates at which to calculate the camber line
%   m  = max camber
%   p  = location of max camber over chord
%   t  = maximum thickness over chord
%   c0 = chord length parameter used to generate the geometry via NACA equations
%   c  = actual length of the undeformed wing; different than c0 because, for example,
%           c0=1.0 results in top and bottom airfoil surfaces meeting at c=1.00893
%   a  = angle of attack
%%%

    % Uniformly-spaced x-coordinates.
    xu = linspace(0, c, 25);

    % Change sign of a so it specifies a clockwise rotation of the airfoil, convert to
    % radians.
    a = -a * pi / 180;

    yc = nan(1, length(xu));
    
	for i = 1:length(xu)

        % Camber line.
        if xu(i) > p*c0
            yc(i) = (m*(c0-xu(i))/(1-p)^2)*(1-2*p+xu(i)/c0);
        else
            yc(i) = (m*(   xu(i))/(  p)^2)*(  2*p-xu(i)/c0);
        end

        % Rotate camber line coordinate pairs for angle of attack.
        rotmat = [cos(a), -sin(a); sin(a), cos(a)];
        tmp = rotmat * [xu(i); yc(i)];
        xu(i) = tmp(1);
        yc(i) = tmp(2);
        
    end
    
    % Evaluate the camber line along the actual airfoil x-coordinates.
    camber_line = interp1(xu, yc, x);
    
end










