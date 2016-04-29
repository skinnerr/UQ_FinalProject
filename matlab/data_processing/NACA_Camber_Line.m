function [ camber_line ] = NACA_Camber_Line( x, m, p, t, c, a )
%%%
% Calculates the camber line y-coordinates of an airfoil with the given parameters.
%
%   x = coordinates at which to calculate the camber line
%   m = max camber
%   p = location of max camber over chord
%   t = maximum thickness over chord
%   c = chord length
%   a = angle of attack
%%%

    % Change sign of a so it specifies a clockwise rotation of the airfoil, convert to
    % radians.
    a = -a * pi / 180;

    camber_line = nan(1,length(x));
    
	for i = 1:length(x)
        
        % Determine region of point relative to location of max camber.
        xgt = x > p*c;

        % Camber line.
        if xgt
            yc = (m*(c-x(i))/(1-p)^2)*(1-2*p+x(i)/c);
        else
            yc = (m*(  x(i))/(  p)^2)*(  2*p-x(i)/c);
        end

        % Rotate camber line coordinate pairs for angle of attack.
        camber_line(i) = cos(a) * x(i) - sin(a) * x(i);
        
    end
    
end










