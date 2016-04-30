function [] = Plot_Cp()

    Set_Default_Plot_Properties();

    nfiles = 10;
    alpha_threshold = 0;
    [NdatL, NstatL] = Load_Processed_Data('LF-1000', nfiles, alpha_threshold);
    [NdatH, NstatH] = Load_Processed_Data('HF-1000', nfiles, alpha_threshold);

    %%%
    % Compute ID and calculate bi-fiedlity model
    %%%
    
    % Form LF and HF matrices
    U_L = nan(length(NdatL(1).x), length(NdatL));
    for i = 1:length(NdatL)
        U_L(:,i) = NdatL(i).cp;
    end
    U_H = nan(length(NdatH(1).x), length(NdatH));
    for i = 1:length(NdatH)
        U_H(:,i) = NdatH(i).cp;
    end
    
    [m, n] = size(U_L); % Number of rows and cols in the LF data
    r = min(m,n);       % Maximal rank
    tol = 1e-3;         % Approximate tolerance for ID
    
    % generate an approximately low-rank matrix
    [U,S,V] = svd(U_L);
    Sdiag = diag(S(1:r,1:r))/S(1,1);

    % plot the singular values of U_c just to see the rank
    semilogy(Sdiag,'ro')
    title('Singular Values of U_L')

    % Show how matrix ID works
    [P,ix] = matrixID(U_L,tol^2); % P is the coefficient matrix and ix is the basis index
    UL_id = U_L(:,ix) * P;
    UF_id = U_H(:,ix) * P;

    %%%
    % Plot stuff
    %%%
    
    figure();
    hold on;
    X = [NstatL.x_avg_twosurf; flip(NstatL.x_avg_twosurf)];
    Y2 =     [NstatL.cp_avg+3*NstatL.cp_std; ...
         flip(NstatL.cp_avg-3*NstatL.cp_std)];
    fill(X,Y2,[0.95, 1.0, 0.95]);
    Y1 = [NstatL.cp_max; flip(NstatL.cp_min)];
    fill(X,Y1,[0.95, 0.95, 1.0]);
    for i = 1:nfiles
        hdots = plot(NstatL.x_avg_twosurf, NdatL(i).cp, 'k.', 'MarkerSize', 3);
    end
    hrange = plot(NstatL.x_avg_twosurf, NstatL.cp_max, 'b-');
             plot(NstatL.x_avg_twosurf, NstatL.cp_min, 'b-');
    h3std  = plot(NstatL.x_avg_twosurf, NstatL.cp_avg+3*NstatL.cp_std, 'g-');
             plot(NstatL.x_avg_twosurf, NstatL.cp_avg-3*NstatL.cp_std, 'g-');
    hmean  = plot(NstatL.x_avg_twosurf, NstatL.cp_avg, 'r-');
    ylim([-3,1.5]);
    xlabel('Streamwise Location CW from Trailing Edge [1 / Streamwise Chord]');
    ylabel('Coefficient of Pressure');
    legend([hdots, hrange, h3std, hmean], {'Realizations', 'Range', '3 stdev', 'Mean'});

end







