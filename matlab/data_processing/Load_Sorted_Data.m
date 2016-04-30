function [ NACA_data ] = Load_Sorted_Data( identifier_string, nfiles, alpha_threshold )

    fprintf('Reading %i files from %s\n',nfiles,identifier_string);

    if strcmp(identifier_string, 'LF-1000') || strcmp(identifier_string, 'HF-1000')
        starting_run_number = 1000;
        params_file = ['/home/ryan/Documents/CU_Boulder/Year_2_2015-2016/', ...
                       'ASEN_6519_Uncertainty_Quantification/Final_Project/', ...
                       'matlab/naca_params_posalpha.csv'];
        if strcmp(identifier_string, 'LF-1000')
            path = ['/home/ryan/Documents/CU_Boulder/Year_2_2015-2016/', ...
                    'ASEN_6519_Uncertainty_Quantification/Final_Project/', ...
                    'data/low-fidelity/'];
        elseif strcmp(identifier_string, 'HF-1000')
            path = ['/home/ryan/Documents/CU_Boulder/Year_2_2015-2016/', ...
                    'ASEN_6519_Uncertainty_Quantification/Final_Project/', ...
                    'data/high-fidelity/'];     
        end 
    end

    NACA_params = csvread(params_file,1,0);
    NACA_data = Load_NACA_Directory(path, 'wing', nfiles);
    NACA_data = Sort_NACA_Data(NACA_params, NACA_data, starting_run_number, ...
                               alpha_threshold);
	if alpha_threshold <= 0
        fprintf('No files removed due to alpha constraint, because threshold <= 0\n');
    else
        fprintf('After removing alpha > %.2f, %i files remain\n',alpha_threshold,nfiles);
    end

end

