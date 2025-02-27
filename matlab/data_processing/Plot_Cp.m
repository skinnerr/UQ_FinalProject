function [] = Plot_Cp()

    Set_Default_Plot_Properties();

    case_id = 5;
    reload_data = false; %true;
    
    if reload_data
        nfiles = 300;
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
            case 4
                [NdatL, NstatL] = Load_Processed_Data('geometryattack-LF', nfiles, alpha_threshold);
                    save('geometryattack-LF-dat.mat', 'NdatL');
                    save('geometryattack-LF-stat.mat', 'NstatL');
                [NdatH, NstatH] = Load_Processed_Data('geometryattack-HF', nfiles, alpha_threshold);
                    save('geometryattack-HF-dat.mat', 'NdatH');
                    save('geometryattack-HF-stat.mat', 'NstatH');
            case 5
                [NdatL, NstatL] = Load_Processed_Data('CW-166k', nfiles, alpha_threshold);
                    save('CW-166k-dat.mat', 'NdatL');
                    save('CW-166k-stat.mat', 'NstatL');
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
            case 4
                load geometryattack-LF-dat.mat
                load geometryattack-LF-stat.mat
                load geometryattack-HF-dat.mat
                load geometryattack-HF-stat.mat
            case 5
                load CW-166k-dat.mat
                load CW-166k-stat.mat
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
    
    second_data_exists = false;
    if second_data_exists
    
        UH = nan(length(NdatH(1).x), length(NdatH));
        for i = 1:length(NdatH)
            UH(:,i) = NdatH(i).cp;
        end

        [m, n] = size(UL); % Number of rows and cols in the LF data
        r = min(m,n);      % Maximal rank
        tol = 1e-1;        % Approximate tolerance for ID

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

        % Compute reduced-order models
        UL_id = UL(:,ix) * P;
        UH_id = UH(:,ix) * P;

        % Compute HF data interpolated at LF coordinates.
        UH_interp = nan(size(UL));
        for i = 1:size(UH,2)
            % Interpolate upper and lower surfaces separately (otherwise interp1 breaks).
            [~,indL] = max(NdatL(i).x);
            [~,indH] = max(NdatH(i).x);
            UH_interp(     1:indL,i) = interp1_bounded(NdatH(i).x(1:indH), ...
                                                       NdatH(i).cp(1:indH), ...
                                                       NdatL(i).x(1:indL));
            UH_interp(indL+1:end ,i) = interp1_bounded(NdatH(i).x(indH+1:end), ...
                                                       NdatH(i).cp(indH+1:end), ...
                                                       NdatL(i).x(indL+1:end));
    %         figure();
    %         hold on;
    %         plot(NdatH(i).x, NdatH(i).cp, 'ko');
    %         plot(xL, UH_interp(:,i), 'r+');
    %         return
        end

        % Compute errors.
        err_id_L = norm(UL - UL_id,'fro')/norm(UL,'fro');
        err_id_H = norm(UH - UH_id,'fro')/norm(UH,'fro');
        err_LvH  = norm(UL - UH_interp,'fro')/norm(UH_interp,'fro');

        fprintf('*** RESULTS ***\n')
        fprintf(' Desired accuracy                = %e\n', tol)
        fprintf(' Accuracy of LF model vs HF      = %e\n', err_LvH)
        fprintf(' Accuracy of ID for coarse model = %e\n', err_id_L)
        fprintf(' Accuracy of ID for fine model   = %e\n', err_id_H)
        fprintf(' Approximation rank              = %d\n', length(ix));

        % Truncate matrices to be plotted if desired.
        if false
            trunc_number = 10;
            UL = UL(:,1:trunc_number);
            UH = UH(:,1:trunc_number);
            UL_id = UL_id(:,1:trunc_number);
            UH_id = UH_id(:,1:trunc_number);
            UH_interp = UH_interp(:,1:trunc_number);
        end
        
    end

    %%%
    % Plot results
    %%%
    
    % Shorthand for average positions.
    xLavg = NstatL.x_avg_twosurf;
    if second_data_exists
    xHavg = NstatH.x_avg_twosurf;
    end
    
    % Shorthand for individual positions.
    xL = nan(size(UL));
    for i = 1:size(UL,2)
        xL(:,i) = NdatL(i).x_twosurf;
    end
    xL_exact = nan(size(UL));
    for i = 1:size(UL,2)
        xL_exact(:,i) = NdatL(i).xnorm;
    end
    if second_data_exists
    xH = nan(size(UH));
    for i = 1:size(UH,2)
        xH(:,i) = NdatH(i).x_twosurf;
    end
    end
    
    rows = 1;
    cols = 2;
    yrange = [min(min(UL)),max(max(UL))];
    yrange(1) = yrange(1) - 0.05 * (yrange(2)-yrange(1));
    yrange(2) = yrange(2) + 0.05 * (yrange(2)-yrange(1));
    figure();
    
    if second_data_exists
    subplot(rows, cols, 1);
    end
    hold on;
    ha = plot(xL, UL, 'ro');
    if second_data_exists
    hb = plot(xL, UL_id, 'k.');
         plot(xL, UH_interp, 'g+');
    hc = plot(xLavg, NstatL.cp_avg, 'c-',  'LineWidth', 2);
    hd = plot(xLavg, mean(UL_id,2), 'm--', 'LineWidth', 2);
    he = plot(xLavg, 5*NstatL.cp_var, 'b-',  'LineWidth', 2);
    hf = plot(xLavg, 5*var(UL_id,0,2),'g--', 'LineWidth', 2);
    hleg = legend([ha(1),hc,he,hb(1),hd,hf], {'Original LF Data', ...
                                              '    mean', ...
                                              '    var*5', ...
                                              sprintf('ID (relerr %.2e)',err_id_L), ...
                                              '    mean', ...
                                              '    var*5'});
    set(hleg, 'Location', 'southeast');
    end
    ylim(yrange);
    ylabel('Cp');
    xlabel('Location on Airfoil');
    set(gca,'YDir','reverse');
    if second_data_exists
    title(sprintf('Case %1i: LF, rank = %i',case_id,length(ix)));
    end
    
    if second_data_exists
    subplot(rows, cols, 2);
    hold on;
    ha = plot(xH, UH, 'ro');
    hb = plot(xH, UH_id, 'k.');
         plot(xL, UH_interp, 'g+');
    hc = plot(xHavg, NstatH.cp_avg, 'c-',  'LineWidth', 2);
    hd = plot(xHavg, mean(UH_id,2), 'm--', 'LineWidth', 2);
    he = plot(xHavg, 5*NstatH.cp_var, 'b-',  'LineWidth', 2);
    hf = plot(xHavg, 5*var(UH_id,0,2),'g--', 'LineWidth', 2);
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
    set(gca,'YDir','reverse');
    title(sprintf('Case %1i: HF, rank = %i',case_id,length(ix)));
    end
    
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

    % Low-fidelity standard Cp plot
    figure();
    hold on;
%     errorbar(NstatL.x_avg, NstatL.cp_avg, ...
%              abs(NstatL.cp_avg - NstatL.cp_min), ...
%              abs(NstatL.cp_avg - NstatL.cp_max));
    ha = plot(xL_exact, UL, 'ro');
%              NstatL.cp_std);
    title('Low-fidelity');
    xlabel('x [m]');
    ylabel('Cp');
    xlim([-0.05,1.05]);
    ylim(yrange);
    set(gca,'YDir','reverse');
    
    if second_data_exists
    % High-fidelity standard Cp plot
    figure();
    hold on;
    errorbar(NstatH.x_avg, NstatH.cp_avg, ...
             abs(NstatH.cp_avg - NstatH.cp_min), ...
             abs(NstatH.cp_avg - NstatH.cp_max));
%              NstatH.cp_std);
    title('High-fidelity');
    xlabel('x [m]');
    ylabel('Cp');
    xlim([-0.05,1.05]);
    ylim(yrange);
    set(gca,'YDir','reverse');
         
    % Bi-fidelity standard Cp plot
    figure();
    hold on;
    errorbar(NstatH.x_avg, mean(UH_id,2), ...
             abs(mean(UH_id,2) - max(UH_id,[],2)), ...
             abs(mean(UH_id,2) - min(UH_id,[],2)));
%              std(UH_id,1,2));
    title('Bi-fidelity');
    xlabel('x [m]');
    ylabel('Cp');
    xlim([-0.05,1.05]);
    ylim(yrange);
    set(gca,'YDir','reverse');
    end

end







