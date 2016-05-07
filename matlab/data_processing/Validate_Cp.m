function [] = Validate_Cp()

    Set_Default_Plot_Properties();

    NACA_params = [0.04, 0.4, 0.12, 1.0, 0.0];
    path = ['/home/ryan/Documents/CU_Boulder/Year_2_2015-2016/', ...
            'ASEN_6519_Uncertainty_Quantification/Final_Project/data/base-case-hf'];
    dat = Load_NACA_Directory(path, 'wing', 0);
    dat = Sort_NACA_Data(NACA_params, dat, 0, 0);
    
    trailing_x = find(dat.x == max(dat.x),1);
    x_special = [dat.x(1:trailing_x); max(dat.x)*2 - dat.x(trailing_x+1:end)];
    x_special = x_special - max(dat.x);

    %%%
    % Plot results
    %%%
    
    figure();
    hold on;
%     hdots = plot(dat.xnorm, dat.y, 'k-o', 'MarkerSize', 3);
    hdots = plot(dat.x, dat.cp, 'r-x', 'MarkerSize', 3);
%     hdots = plot(x_special, dat.cp, 'r-x', 'MarkerSize', 3);
    xlabel('Streamwise Location CW from Trailing Edge [1 / Streamwise Chord]');
    ylabel('Coefficient of Pressure');
    legend([hdots], {'CFD'});

end







