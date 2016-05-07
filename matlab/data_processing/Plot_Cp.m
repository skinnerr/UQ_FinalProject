function [] = Plot_Cp()

    Set_Default_Plot_Properties();

    case_id = 1; % Valid numbers: 1,2,3
    reload_data = false;
    
    if reload_data
        nfiles = 1000;
        alpha_threshold = 0;
        switch case_id
            case 1
                [NdatL, NstatL] = Load_Processed_Data('original-LF', nfiles, alpha_threshold);
                    save('original-LF-dat.mat', 'NdatL');
                    save('original-LF-stat.mat', 'NstatL');
                [NdatH, NstatH] = Load_Processed_Data('original-HF', nfiles, alpha_threshold);
                    save('original-HF-dat.mat', 'NdatH');
                    save('original-HF-stat.mat', 'NstatH');
            case 2
                [NdatL, NstatL] = Load_Processed_Data('geometry-LF', nfiles, alpha_threshold);
                    save('geometry-LF-dat.mat', 'NdatL');
                    save('geometry-LF-stat.mat', 'NstatL');
                [NdatH, NstatH] = Load_Processed_Data('geometry-HF', nfiles, alpha_threshold);
                    save('geometry-HF-dat.mat', 'NdatH');
                    save('geometry-HF-stat.mat', 'NstatH');
            case 3
                [NdatL, NstatL] = Load_Processed_Data('attack-LF', nfiles, alpha_threshold);
                    save('attack-LF-dat.mat', 'NdatL');
                    save('attack-LF-stat.mat', 'NstatL');
                [NdatH, NstatH] = Load_Processed_Data('attack-HF', nfiles, alpha_threshold);
                    save('attack-HF-dat.mat', 'NdatH');
                    save('attack-HF-stat.mat', 'NstatH');
        end
    else
        switch case_id
            case 1
                load original-LF-dat.mat
                load original-LF-stat.mat
                load original-HF-dat.mat
                load original-HF-stat.mat
            case 2
                load geometry-LF-dat.mat
                load geometry-LF-stat.mat
                load geometry-HF-dat.mat
                load geometry-HF-stat.mat
            case 3
                load attack-LF-dat.mat
                load attack-LF-stat.mat
                load attack-HF-dat.mat
                load attack-HF-stat.mat
        end
    end

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
    tol = 1e-10;         % Approximate tolerance for ID
    
    % generate an approximately low-rank matrix
    [~,SL,~] = svd(UL);
    SLdiag = diag(SL(1:r,1:r))/SL(1,1);
    [~,SH,~] = svd(UH);
    SHdiag = diag(SH(1:r,1:r))/SH(1,1);

    % plot the singular values of U_c just to see the rank
    figure();
    
    subplot(1,2,1);
    semilogy(SLdiag,'ro');
    axis tight;
    tmp = get(gca, 'YLim');
    xlabel('Index');
    ylabel('Singular Value Magnitude');
    title(sprintf('Case %1i: LF',case_id));
    
    subplot(1,2,2);
    semilogy(SHdiag,'ro');
    axis tight;
    ylim(tmp);
    xlabel('Index');
    title(sprintf('Case %1i: HF',case_id));

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
    yrange = [-3.5,1.5];
    figure();
    
    subplot(rows, cols, 1);
    hold on;
    ha = plot(xL, UL, 'ro');
    hb = plot(xL, UL_id, 'k.');
    hc = plot(xL, NstatL.cp_avg, 'c-',  'LineWidth', 2);
    hd = plot(xL, mean(UL_id,2), 'm--', 'LineWidth', 2);
    he = plot(xL, 5*NstatL.cp_var, 'b-',  'LineWidth', 2);
    hf = plot(xL, 5*var(UL_id,0,2),'g--', 'LineWidth', 2);
    hleg = legend([ha(1),hc,he,hb(1),hd,hf], {'Original LF Data', ...
                                              '    mean', ...
                                              '    var*5', ...
                                              sprintf('ID (relerr %.2e)',err_id_L), ...
                                              '    mean', ...
                                              '    var*5'});
    set(hleg, 'Location', 'southeast');
    ylim(yrange);
    ylabel('Cp');
    xlabel('Location on Airfoil');
    title(sprintf('Case %1i: LF, rank = %i',case_id,length(ix)));
    
    subplot(rows, cols, 2);
    hold on;
    ha = plot(xH, UH, 'ro');
    hb = plot(xH, UH_id, 'k.');
    hc = plot(xH, NstatH.cp_avg, 'c-',  'LineWidth', 2);
    hd = plot(xH, mean(UH_id,2), 'm--', 'LineWidth', 2);
    he = plot(xH, 5*NstatH.cp_var, 'b-',  'LineWidth', 2);
    hf = plot(xH, 5*var(UH_id,0,2),'g--', 'LineWidth', 2);
    hleg = legend([ha(1),hc,he,hb(1),hd,hf], {'Original HF Data', ...
                                              '    mean', ...
                                              '    var*5', ...
                                              sprintf('Bi-Fidelity (relerr %.2e)',err_id_H), ...
                                              '    mean', ...
                                              '    var*5'});
    set(hleg, 'Location', 'southeast');
    ylim(yrange);
    ylabel('Cp');
    xlabel('Location on Airfoil');
    title(sprintf('Case %1i: HF, rank = %i',case_id,length(ix)));
    
%     figure();
%     hold on;
%     X = [NstatL.x_avg_twosurf; flip(NstatL.x_avg_twosurf)];
%     Y2 =     [NstatL.cp_avg+3*NstatL.cp_std; ...
%          flip(NstatL.cp_avg-3*NstatL.cp_std)];
%     fill(X,Y2,[0.95, 1.0, 0.95]);
%     Y1 = [NstatL.cp_max; flip(NstatL.cp_min)];
%     fill(X,Y1,[0.95, 0.95, 1.0]);
%     for i = 1:length(NdatL)
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







