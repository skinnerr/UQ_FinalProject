function [] = Plot_Cp()

    nfiles = 1;

    %%%
    % Load data from NACA realizations
    %%%
    
    path = ['/home/ryan/Documents/CU_Boulder/Year_2_2015-2016/', ...
            'ASEN_6519_Uncertainty_Quantification/Final_Project/data/high-fidelity/'];
             
    NACA_data = Load_NACA_Directory(path, 'wing', nfiles);
    
    %%%
    % Sort points going around airfoil circumference
    %%%
    
    params_file = ['/home/ryan/Documents/CU_Boulder/Year_2_2015-2016/', ...
                   'ASEN_6519_Uncertainty_Quantification/Final_Project/', ...
                   'matlab/naca_params_posalpha.csv'];
    
	starting_run_number = 1000;
    NACA_data = Sort_NACA_Data(params_file, NACA_data, starting_run_number);

    %%%
    % Plot stuff!
    %%%
    
    figure();
    for i = 1:nfiles
        hold on;
        plot(NACA_data(i).xnorm, NACA_data(i).y, '-');
    end

end