function [] = Plot_Cp()

    nfiles = 100;

    %%%
    % Load data from NACA realizations
    %%%
    
    path = ['/home/ryan/Documents/CU_Boulder/Year_2_2015-2016/', ...
                 'ASEN_6519_Uncertainty_Quantification/Final_Project/data/high-fidelity/'];
    NACA_data = Load_NACA_Directory(path, 'wing', nfiles);
    
    %%%
    % Sort points going around airfoil circumference
    %%%
    
    NACA_data = Sort_NACA_Data(path, NACA_data);

    %%%
    % Plot stuff!
    %%%
    
    figure();
    for i = 1:nfiles
        hold on;
        plot(NACA_data(i).xnorm, NACA_data(i).cp, '.');
    end

end