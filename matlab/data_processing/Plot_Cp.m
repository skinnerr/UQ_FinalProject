function [] = Plot_Cp()

    Set_Default_Plot_Properties();

    nfiles = 100;
    alpha_threshold = 0;
    NACA_data = Load_Sorted_Data('LF-1000', nfiles, alpha_threshold);

    %%%
    % Plot stuff!
    %%%
    
    x_avg  = zeros(length(NACA_data(1).xnorm),1);
    cp_max = NACA_data(1).cp;
    cp_min = cp_max;
    cp_avg = zeros(length(cp_max),1);
    cp_var = zeros(length(NACA_data(1).xnorm),1);
    for i = 1:nfiles
        x_avg = x_avg + NACA_data(i).xnorm / nfiles;
        cp_avg = cp_avg + NACA_data(i).cp / nfiles;
        cp_var = cp_var + (NACA_data(i).cp).^2 / nfiles;
        ind = NACA_data(i).cp > cp_max;
        cp_max(ind) = NACA_data(i).cp(ind);
        ind = NACA_data(i).cp < cp_min;
        cp_min(ind) = NACA_data(i).cp(ind);
    end
    cp_var = cp_var - cp_avg.^2;
    cp_std = sqrt(cp_var);
    
    % Create x_avg_twosurf coordinate, which has x<0 on upper surface, x>0 on lower
    % surface, and x=0 on the trailing edge.
    trailing_x = find(x_avg == max(x_avg),1);
    x_avg_twosurf = [x_avg(1:trailing_x); max(x_avg)*2 - x_avg(trailing_x+1:end)];
    x_avg_twosurf = x_avg_twosurf - max(x_avg);
    
    figure();
    hold on;
    X = [x_avg_twosurf; flip(x_avg_twosurf)];
    Y2 = [cp_avg+3*cp_std; flip(cp_avg-3*cp_std)];
    fill(X,Y2,[0.95, 1.0, 0.95]);
    Y1 = [cp_max; flip(cp_min)];
    fill(X,Y1,[0.95, 0.95, 1.0]);
    for i = 1:nfiles
        hdots = plot(x_avg_twosurf, NACA_data(i).cp, 'k.', 'MarkerSize', 3);
    end
    hrange = plot(x_avg_twosurf, cp_max, 'b-');
             plot(x_avg_twosurf, cp_min, 'b-');
    h3std  = plot(x_avg_twosurf, cp_avg+3*cp_std, 'g-');
             plot(x_avg_twosurf, cp_avg-3*cp_std, 'g-');
    hmean  = plot(x_avg_twosurf, cp_avg, 'r-');
    ylim([-3,1.5]);
    xlabel('Streamwise Location CW from Trailing Edge [1 / Streamwise Chord]');
    ylabel('Coefficient of Pressure');
    legend([hdots, hrange, h3std, hmean], {'Realizations', 'Range', '3 stdev', 'Mean'});

end







