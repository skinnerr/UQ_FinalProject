function [] = Plot_Cp()

    Set_Default_Plot_Properties();

    nfiles = 100;
    alpha_threshold = 0;
    [Ndat, Nstat] = Load_Processed_Data('LF-1000', nfiles, alpha_threshold);

    %%%
    % Plot stuff!
    %%%
    
    figure();
    hold on;
    X = [Nstat.x_avg_twosurf; flip(Nstat.x_avg_twosurf)];
    Y2 =     [Nstat.cp_avg+3*Nstat.cp_std; ...
         flip(Nstat.cp_avg-3*Nstat.cp_std)];
    fill(X,Y2,[0.95, 1.0, 0.95]);
    Y1 = [Nstat.cp_max; flip(Nstat.cp_min)];
    fill(X,Y1,[0.95, 0.95, 1.0]);
    for i = 1:nfiles
        hdots = plot(Nstat.x_avg_twosurf, Ndat(i).cp, 'k.', 'MarkerSize', 3);
    end
    hrange = plot(Nstat.x_avg_twosurf, Nstat.cp_max, 'b-');
             plot(Nstat.x_avg_twosurf, Nstat.cp_min, 'b-');
    h3std  = plot(Nstat.x_avg_twosurf, Nstat.cp_avg+3*Nstat.cp_std, 'g-');
             plot(Nstat.x_avg_twosurf, Nstat.cp_avg-3*Nstat.cp_std, 'g-');
    hmean  = plot(Nstat.x_avg_twosurf, Nstat.cp_avg, 'r-');
    ylim([-3,1.5]);
    xlabel('Streamwise Location CW from Trailing Edge [1 / Streamwise Chord]');
    ylabel('Coefficient of Pressure');
    legend([hdots, hrange, h3std, hmean], {'Realizations', 'Range', '3 stdev', 'Mean'});

end







