function [] = Plot_Cp()

    Set_Default_Plot_Properties();

    nfiles = 100;
    alpha_threshold = 0;
    [NdatL, NstatL] = Load_Processed_Data('original-LF', nfiles, alpha_threshold);
    [NdatH, NstatH] = Load_Processed_Data('original-HF', nfiles, alpha_threshold);
%     [NdatL, NstatL] = Load_Processed_Data('geometry-LF', nfiles, alpha_threshold);
%     [NdatH, NstatH] = Load_Processed_Data('geometry-HF', nfiles, alpha_threshold);
%     [NdatL, NstatL] = Load_Processed_Data('attack-LF', nfiles, alpha_threshold);
%     [NdatH, NstatH] = Load_Processed_Data('attack-HF', nfiles, alpha_threshold);

    %%%
    % Compute ID and calculate bi-fiedlity model
    %%%
    
    % Form LF and HF matrices
    UL = nan(length(NdatL(1).x), length(NdatL));
    for i = 1:length(NdatL)
        UL(:,i) = NdatL(i).cp;
    end
    UH = nan(length(NdatH(1).x), length(NdatH));
    for i = 1:length(NdatH)
        UH(:,i) = NdatH(i).cp;
    end
    
    [m, n] = size(UL); % Number of rows and cols in the LF data
    r = min(m,n);       % Maximal rank
    tol = 0.8;         % Approximate tolerance for ID
    
    % generate an approximately low-rank matrix
    [~,SL,~] = svd(UL);
    SLdiag = diag(SL(1:r,1:r))/SL(1,1);
    [~,SH,~] = svd(UH);
    SHdiag = diag(SH(1:r,1:r))/SH(1,1);

    % plot the singular values of U_c just to see the rank
    figure();
    subplot(1,2,1);
    semilogy(SLdiag,'ro')
    title('Singular Values of U_L')
    subplot(1,2,2);
    semilogy(SHdiag,'ro')
    title('Singular Values of U_H')

    % Perform ID and bi-fidelity modeling
    [P,ix] = matrixID(UL,tol^2); % P is the coefficient matrix and ix is the basis index
    
    % Truncate the expansion and compute reduced-order models
    UL_id = UL(:,ix) * P;
    UH_id = UH(:,ix) * P;
    
    % Compute errors.
    err_id_L = norm(UL - UL_id,'fro')/norm(UL,'fro');
    err_id_H = norm(UH - UH_id,'fro')/norm(UH,'fro');

    fprintf('*** RESULTS ***\n')
    fprintf(' Desired accuracy                = %e\n', tol)
    fprintf(' Accuracy of ID for coarse model = %e\n', err_id_L) 
    fprintf(' Accuracy of ID for fine model   = %e\n', err_id_H) 
    fprintf(' Approximation rank              = %d\n', length(ix));

    %%%
    % Plot results
    %%%
    
    xL = NstatL.x_avg_twosurf;
    xH = NstatH.x_avg_twosurf;
    
    rows = 1;
    cols = 2;
    yrange = [-2,1];
    figure();
    
    subplot(rows, cols, 1);
    hold on;
    ha = plot(xL, UL, 'ro');
    hb = plot(xL, UL_id, 'k.');
    hleg = legend([ha(1),hb(1)], {'Original LF Data','ID Approximation of LF'});
    set(hleg, 'Location', 'southeast');
    ylim(yrange);
    ylabel('Cp');
    xlabel('Streamwise Location');
    title('Low-Fidelity');
    
    subplot(rows, cols, 2);
    hold on;
    ha = plot(xH, UH, 'ro');
    hb = plot(xH, UH_id, 'k.');
    hleg = legend([ha(1),hb(1)], {'Original HF Data','Bi-Fidelity Model'});
    set(hleg, 'Location', 'southeast');
    ylim(yrange);
    ylabel('Cp');
    xlabel('Streamwise Location');
    title('High-Fidelity');
    
%     figure();
%     hold on;
%     X = [NstatL.x_avg_twosurf; flip(NstatL.x_avg_twosurf)];
%     Y2 =     [NstatL.cp_avg+3*NstatL.cp_std; ...
%          flip(NstatL.cp_avg-3*NstatL.cp_std)];
%     fill(X,Y2,[0.95, 1.0, 0.95]);
%     Y1 = [NstatL.cp_max; flip(NstatL.cp_min)];
%     fill(X,Y1,[0.95, 0.95, 1.0]);
%     for i = 1:nfiles
%         hdots = plot(NstatL.x_avg_twosurf, NdatL(i).cp, 'k.', 'MarkerSize', 3);
%     end
%     hrange = plot(NstatL.x_avg_twosurf, NstatL.cp_max, 'b-');
%              plot(NstatL.x_avg_twosurf, NstatL.cp_min, 'b-');
%     h3std  = plot(NstatL.x_avg_twosurf, NstatL.cp_avg+3*NstatL.cp_std, 'g-');
%              plot(NstatL.x_avg_twosurf, NstatL.cp_avg-3*NstatL.cp_std, 'g-');
%     hmean  = plot(NstatL.x_avg_twosurf, NstatL.cp_avg, 'r-');
%     ylim([-3,1.5]);
%     xlabel('Streamwise Location CW from Trailing Edge [1 / Streamwise Chord]');
%     ylabel('Coefficient of Pressure');
%     legend([hdots, hrange, h3std, hmean], {'Realizations', 'Range', '3 stdev', 'Mean'});

end







