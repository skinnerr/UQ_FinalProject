function [ ] = naca_realizations()
    
    Set_Default_Plot_Properties();
    
    % For a NACA airfoil numbered xyzz...
    m = getRange(0.04, 0); % Max camber (x/100)
    p = getRange(0.4, 0); % Location of max camber (y/10)
    t = getRange(0.12, 0); % Max thickness over chord (zz / 100)
    a = [-2,2];  % Angle of attack
    
    points = 201;
    c = 1.0;        % Chord length
    x = linspace(0, c, points);
    
    % Plot airfoil realizations
    figure();
    hold on;
    
    N = 1000;
    params = nan(N,5);
    for i = 1:N
        rm = randu(m);
        rp = randu(p);
        rt = randu(t);
        rc = c;
        ra = randu(a);
        [yt, xu, yu, xl, yl] = naca_coords(x, rm, rp, rt, ra);
        params(i,:) = [rm, rp, rt, rc, ra];
        if i < 10
            plot([xu,flip(xl)], [yu,flip(yl)]);
        end
    end
    plot( x,  yt, ':k', 'LineWidth', 3);
    plot( x, -yt, ':k', 'LineWidth', 3);
    axis square;
    axis equal;
    
    writetable(array2table(params, 'VariableNames', {'m','p','t','c','a'}), ...
               '../naca_params_4412alphaneg2to2.csv');
    
    
    if false
        % Plot the un-cambered airfoil.
        figure();
        plot( [x,flip(x)], [yt,flip(-yt)]);
        axis equal;
        xlabel('x');
        ylabel('y');
        xlim([0,c]+[-1,1]*c/10);
        ylim(1.1*[-t(2),t(2)]*c);
    end
    
end

function [ x ] = randu( range )

    dif = range(2) - range(1);
    x = range(1) + rand() * dif;

end

function [ range ] = getRange( val, pmerr )

    range = [val - val*pmerr, val + val*pmerr];

end