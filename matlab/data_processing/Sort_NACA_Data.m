function [ data_out ] = Sort_NACA_Data( naca_params, data, starting_run_number, thresh )
%%%
% Sorts data points within each entry of the NACA data structure, such that they progress
% around the circumference of the deformed wing. This is done by computing the camber line
% of the deformed airfoil, and identifying nodes above and below it. These two sets of
% nodes are then sorted by increasing and decreasing x-coordinate, respectively.
%%%
    
    try
        
        load_str = 'Sorting nodal point data...';
        hwait = waitbar(0,sprintf('%s %i / %i',load_str,1,length(data)), ...
                        'CreateCancelBtn','setappdata(gcbf,''cancel_loading'',1)');
        setappdata(hwait,'cancel_loading',0);
        
        n_skipped = 0;
        for i = 1:length(data)
            
            if getappdata(hwait,'cancel_loading')
                delete(hwait);
                error('Nodal point sorting canceled by user. Program stopped.');
            end
            waitbar(i/length(data),hwait,sprintf('%s %i / %i',load_str,i,length(data)));
            
            % Extract the run number from the loaded file's ID string (from its filename).
            run_number = strsplit(data(i).IDstr,'-');
            run_number = str2num(run_number{end}); %#ok<ST2NM>
            run_number = run_number - starting_run_number;
            
            % Grab the parameters used to generate this realization of airfoil data.
            m  = naca_params(run_number,1);
            p  = naca_params(run_number,2);
            c0 = naca_params(run_number,4); % Value of c used to generate the geometry.
            c  = 1.00893; % Actual length of the geometry.
            a  = naca_params(run_number,5);
            
            % Skip if angle of attack is above a certain threshold.
            if a > thresh
                n_skipped = n_skipped + 1;
                continue;
            end
            
            % Sort x-coordinates.
            [x_sorted, order] = sort(data(i).x);
            
            % Calculate the camber line for this airfoil.
            camber_line = NACA_Camber_Line(x_sorted, m, p, c0, c, a);
            
            % Determine what points are on top or bottom surface of airfoil.
            top_indices = data(i).y(order) > camber_line;
            bot_indices = ~top_indices;
            
            % Split the order to top indices go first, followed by negative indies in
            % reverse.
            order = [order(top_indices); flip(order(bot_indices))];
            
            % Re-order the elements of the data matrix corresponding to this run.
            data_out(i-n_skipped).IDstr = data(i).IDstr; %#ok<*AGROW>
            data_out(i-n_skipped).x     = data(i).x(order);
            data_out(i-n_skipped).xnorm = data(i).xnorm(order);
            data_out(i-n_skipped).y     = data(i).y(order);
            data_out(i-n_skipped).p     = data(i).p(order);
            data_out(i-n_skipped).cp    = data(i).cp(order);
            data_out(i-n_skipped).ref   = data(i).ref;
            data_out(i-n_skipped).filename = data(i).filename;
            
        end
        
    catch err
        
        delete(hwait);
        rethrow(err);
        
    end

    delete(hwait);

end

    
    
    
    
    