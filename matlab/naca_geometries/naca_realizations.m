function [ ] = naca_realizations()
    
    Set_Default_Plot_Properties();
    
    % For a NACA airfoil numbered xyzz...
    m = [0.0, 0.1]; % Max camber (x/100)
    p = [0.3, 0.6]; % Location of max camber (y/10)
    t = [0.05, 0.2]; % Max thickness over chord (zz / 100)
    a = [0, 7];  % Angle of attack
    
    points = 201;
    c = 1.0;        % Chord length
    x = linspace(0, c, points);
    
    % Plot airfoil realizations
    figure();
    hold on;
    
    N = 10000;
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
    
    writetable(array2table(params, 'VariableNames', {'m','p','t','c','a'}), 'naca_params_posalpha2.csv');
    
    
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