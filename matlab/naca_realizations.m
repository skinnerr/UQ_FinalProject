function [ ] = naca_realizations()
    
    Set_Default_Plot_Properties();
    
    % For a NACA airfoil numbered xyzz...
    m = 4/100;  % Max camber (x/100)
    p = 4/10;   % Location of max camber (y/10)
    t = 12/100; % Max thickness over chord (zz / 100)
    c = 1.0;    % Chord length
    points = 201;
    x = linspace(0, c, points);
    
    % Plot airfoil realizations
    figure();
    hold on;
    
    N = 1;
    for i = 1:N
        [yt, xu, yu, xl, yl] = naca_coords(x, m*(1+randn), p, t);
        plot([xu,flip(xl)], [yu,flip(yl)]);
    end
    plot( x, yt, ':k');
    axis square;
    axis equal;
    
end