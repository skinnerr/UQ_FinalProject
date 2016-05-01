function [] = Validate_Cp()

    Set_Default_Plot_Properties();

    NACA_params = [0.04, 0.4, 0.12, 1.0, 0.0];
    path = ['/home/ryan/Documents/CU_Boulder/Year_2_2015-2016/', ...
            'ASEN_6519_Uncertainty_Quantification/Final_Project/data/base-case-hf'];
    dat = Load_NACA_Directory(path, 'wing', 0);
    dat = Sort_NACA_Data(NACA_params, dat, 0, 0);

    %%%
    % Plot results
    %%%
    
    figure();
    hold on;
    hdots = plot(dat.xnorm, dat.cp, 'k-o', 'MarkerSize', 3);
    ylim([-3,1.5]);
    xlabel('Streamwise Location CW from Trailing Edge [1 / Streamwise Chord]');
    ylabel('Coefficient of Pressure');
    legend([hdots], {'CFD'});

end







