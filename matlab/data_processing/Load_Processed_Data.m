function [ NACA_data, NACA_stats ] = Load_Processed_Data( identifier_string, nfiles, ...
                                                          alpha_threshold )

    fprintf('Reading %i files from %s\n',nfiles,identifier_string);
    
    main_data_dir = ['/home/ryan/Documents/CU_Boulder/Year_2_2015-2016/', ...
                    'ASEN_6519_Uncertainty_Quantification/Final_Project/data/'];

    if strcmp(identifier_string, 'original-LF') || strcmp(identifier_string, 'original-HF')
        starting_run_number = 1000;
        params_file = [main_data_dir, 'naca_params_posalpha.csv'];
        if strcmp(identifier_string, 'original-LF')
            path = [main_data_dir, 'original-lf/'];
        elseif strcmp(identifier_string, 'original-HF')
            path = [main_data_dir, 'original-hf/'];
        end
    elseif strcmp(identifier_string, 'attack-LF') || strcmp(identifier_string, 'attack-HF')
        starting_run_number = 0;
        params_file = [main_data_dir, 'naca_params_4412alphaneg2to2.csv'];
        if strcmp(identifier_string, 'attack-LF')
            path = [main_data_dir, 'attack-lf/'];
        elseif strcmp(identifier_string, 'attack-HF')
            path = [main_data_dir, 'attack-hf/'];
        end
    elseif strcmp(identifier_string, 'geometry-LF') || strcmp(identifier_string, 'geometry-HF')
        starting_run_number = 0;
        params_file = [main_data_dir, 'naca_params_4412at5percentgeomerr.csv'];
        if strcmp(identifier_string, 'geometry-LF')
            path = [main_data_dir, 'geometry-lf/'];
        elseif strcmp(identifier_string, 'geometry-HF')
            path = [main_data_dir, 'geometry-hf/'];
        end
    elseif strcmp(identifier_string, 'geometryattack-LF') || strcmp(identifier_string, 'geometryattack-HF')
        starting_run_number = 0;
        params_file = [main_data_dir, 'naca_params_4412geomerr10alpha2tom2.csv'];
        if strcmp(identifier_string, 'geometryattack-LF')
            path = [main_data_dir, 'geometryattack-lf/'];
        elseif strcmp(identifier_string, 'geometryattack-HF')
            path = [main_data_dir, 'geometryattack-hf/'];
        end
    elseif strcmp(identifier_string, 'CW-166k')
        starting_run_number = 0;
        params_file = [main_data_dir, 'naca_params_4412geomerr15pctalpha13p87.csv'];
        path = [main_data_dir, 'A0a_2face_CW_166k-geomerr15pctalpha13p87/'];
    else
        error('String "%s" not recognized',identifier_string);
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
    
    if nfiles <= 0
        nfiles = length(NACA_data);
    end
    
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
    
    NACA_stats.x_avg = x_avg;
    NACA_stats.cp_max = cp_max;
    NACA_stats.cp_min = cp_min;
    NACA_stats.cp_avg = cp_avg;
    NACA_stats.cp_var = cp_var;
    NACA_stats.cp_std = cp_std;
    NACA_stats.x_avg_twosurf = x_avg_twosurf;

end

