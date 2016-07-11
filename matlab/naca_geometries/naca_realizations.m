function [ ] = naca_realizations()
    
    Set_Default_Plot_Properties();
    
%     % For a NACA airfoil numbered xyzz...
    m0 = 0.04;
    p0 = 0.40;
    t0 = 0.12;
    a0 = 13.87;
    m = getRange(m0, 0.15); % Max camber (x/100)
    p = getRange(p0, 0.15); % Location of max camber (y/10)
    t = getRange(t0, 0.15); % Max thickness over chord (zz / 100)
    a = getRange(a0, 0.00); % Angle of attack
%     a = [-2,2];  % Angle of attack
%     m = getRange(0.04, 0); % Max camber (x/100)
%     p = getRange(0.4, 0); % Location of max camber (y/10)
%     t = getRange(0.12, 0); % Max thickness over chord (zz / 100)
%     a = [-2,2];  % Angle of attack
    
    points = 201;
    c = 1.0;        % Chord length
    x = linspace(0, c, points);
    
    %
    % Generate airfoil realizations
    %
    N = 1000;
    geom_params = nan(N,5);
    geom_coords = nan(N,4,length(x));
    for n = 1:N
        rm = randu(m);
        rp = randu(p);
        rt = randu(t);
        rc = c;
        ra = randu(a);
        [~, xu, yu, xl, yl] = naca_coords(x, rm, rp, rt, ra);
        geom_params(n,:) = [rm, rp, rt, rc, ra];
        geom_coords(n,1,:) = xu;
        geom_coords(n,2,:) = yu;
        geom_coords(n,3,:) = xl;
        geom_coords(n,4,:) = yl;
    end
    
    writetable(array2table(geom_params, 'VariableNames', {'m','p','t','c','a'}), ...
               '../naca_params_4412geomerr15pctalpha13p87.csv');
    
    
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
    
    % Plot ranges for surface coordinates.
    if true
        figure();
        hold on;
        % Calculate distance from test points for upper and lower surfaces.
        du   = nan(N,points);
        dl   = nan(N,points);
        for i = 1:points
        for n = 1:N
            testx = mean(geom_coords(:,1,i));
            testy = -0.1;
            du(n,i) = sqrt( (testx-geom_coords(n,1,i))^2 + (testy-geom_coords(n,2,i))^2 );
            testx = mean(geom_coords(:,3,i));
            testy = 0.2;
            dl(n,i) = sqrt( (testx-geom_coords(n,3,i))^2 + (testy-geom_coords(n,4,i))^2 );
        end
        end
        % Extract coordinates of bounding surfaces.
        % Min/max notation corresponds to distance to test point --- not y-coordinate.
        [~,ind] = min(du);
        minxu = diag(squeeze(geom_coords(ind,1,:)));
        minyu = diag(squeeze(geom_coords(ind,2,:)));
        [~,ind] = max(du);
        maxxu = diag(squeeze(geom_coords(ind,1,:)));
        maxyu = diag(squeeze(geom_coords(ind,2,:)));
        [~,ind] = min(dl);
        minxl = diag(squeeze(geom_coords(ind,3,:)));
        minyl = diag(squeeze(geom_coords(ind,4,:)));
        [~,ind] = max(dl);
        maxxl = diag(squeeze(geom_coords(ind,3,:)));
        maxyl = diag(squeeze(geom_coords(ind,4,:)));
        % Determine where max lower surface and min upper surface cross
        crossind = find(minyl > minyu);
        if crossind
            crossind = crossind(1)-1;
        else
            crossind = points;
        end
        minxl = minxl(1:crossind);
        minyl = minyl(1:crossind);
        minxu = minxu(1:crossind);
        minyu = minyu(1:crossind);
%         plot(minyl,'r');
%         plot(maxyl,'g');
%         plot(minyu,'b');
%         plot(maxyu,'k');
        % Plot upper and lower surface bounds.
        X = [maxxu; flip(maxxl); minxl; flip(minxu)];
        Y = [maxyu; flip(maxyl); minyl; flip(minyu)];
        color = [175 225 225] / 255;
        color = 'r';
        fill(X,Y,color,'EdgeColor','k');
        % Generate and plot NACA 4412 airfoil.
        [~, xu, yu, xl, yl] = naca_coords(x, m0, p0, t0, a0);
        plot(xu,yu,'k');
        plot(xl,yl,'k');
        xlim([min(x),max(x)] + [-1,1]*range(x)/20);
        % Set axes.
        xlabel('x [m]');
        ylabel('y [m]');
        axis square;
        axis equal;
    end
    
end

function [ x ] = randu( range )

    dif = range(2) - range(1);
    x = range(1) + rand() * dif;

end

function [ range ] = getRange( val, pmerr )

    range = [val - val*pmerr, val + val*pmerr];

end