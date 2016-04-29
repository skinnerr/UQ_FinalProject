function [] = Plot_Cp()

    Set_Default_Plot_Properties();

    nfiles = 100;
    fprintf('Reading %i files\n',nfiles);

    %%%
    % Load data from NACA realizations
    %%%
    
    path = ['/home/ryan/Documents/CU_Boulder/Year_2_2015-2016/', ...
            'ASEN_6519_Uncertainty_Quantification/Final_Project/data/low-fidelity/'];
             
    NACA_data = Load_NACA_Directory(path, 'wing', nfiles);
    
    %%%
    % Sort points going around airfoil circumference
    %%%
    
    params_file = ['/home/ryan/Documents/CU_Boulder/Year_2_2015-2016/', ...
                   'ASEN_6519_Uncertainty_Quantification/Final_Project/', ...
                   'matlab/naca_params_posalpha.csv'];
    
	starting_run_number = 1000;
	naca_params = csvread(params_file,1,0);
    alpha_threshold = 3.0;
    NACA_data = Sort_NACA_Data(naca_params, NACA_data, ...
                               starting_run_number, alpha_threshold);
                           
	nfiles = length(NACA_data); % Some could have been filtered out by alpha_threshold.
    fprintf('After removing alpha > %.2f, %i files remain\n',alpha_threshold,nfiles);

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
    
    % Create x_avg_twosurf coordinate, which has x<0 on upper surface, x>0 on lower
    % surface, and x=0 on the trailing edge.
    trailing_x = find(x_avg == max(x_avg),1);
    x_avg_twosurf = [x_avg(1:trailing_x); max(x_avg)*2 - x_avg(trailing_x+1:end)];
    x_avg_twosurf = x_avg_twosurf - max(x_avg);
    
    figure();
    hold on;
    X = [x_avg_twosurf; flip(x_avg_twosurf)];
    Y = [cp_max; flip(cp_min)];
    fill(X,Y,[0.95, 0.95, 1.0]);
    for i = 1:nfiles
        hdots = plot(x_avg_twosurf, NACA_data(i).cp, 'k.', 'MarkerSize', 3);
    end
    hrange = plot(x_avg_twosurf, cp_max, 'b-');
             plot(x_avg_twosurf, cp_min, 'b-');
             plot(x_avg_twosurf, cp_avg+cp_var, 'g-');
             plot(x_avg_twosurf, cp_avg-cp_var, 'g-');
    hmean  = plot(x_avg_twosurf, cp_avg, 'r-');
    ylim([-3,1.5]);
    xlabel('Streamwise Location CW from Trailing Edge [1 / Streamwise Chord]');
    ylabel('Coefficient of Pressure');
    legend([hdots, hrange, hmean], {'Realizations', 'Range', 'Mean'});

end







