function [ data_out ] = Sort_NACA_Data( params_path, data, starting_run_number )
%%%
% Sorts data points within each entry of the NACA data structure, such that they progress
% around the circumference of the deformed wing. This is done by computing the camber line
% of the deformed airfoil, and identifying nodes above and below it. These two sets of
% nodes are then sorted by increasing and decreasing x-coordinate, respectively.
%%%
	
	naca_params = csvread(params_path,1,0);
    
    %%%
    % Load the remaining '.csv' files and build the output data structure.
    %%%
    
%     try
        
        load_str = 'Sorting nodal point data...';
        hwait = waitbar(0,sprintf('%s %i / %i',load_str,1,length(data)), ...
                        'CreateCancelBtn','setappdata(gcbf,''cancel_loading'',1)');
        setappdata(hwait,'cancel_loading',0);
        
        for i = 1:length(data)
            
            if getappdata(hwait,'cancel_loading')
                delete(hwait);
                error('Nodal point sorting canceled by user. Program stopped.');
            end
            waitbar(i/length(data),hwait,sprintf('%s %i / %i',load_str,i,length(data)));
            
            % Extract the run number from the loaded file's ID string (from its filename).
            run_number = strsplit(data(i).IDstr,'-');
            run_number = str2num(run_number{end});
            run_number = run_number - starting_run_number;
            
            % Grab the parameters used to generate this realization of airfoil data.
            m = naca_params(run_number,1);
            p = naca_params(run_number,2);
            t = naca_params(run_number,3);
            c = naca_params(run_number,4);
            a = naca_params(run_number,5);
            c = 1.00893; % Value in the file (1.0) is wrong; actual geometry uses 1.00893.
            
            % Sort x-coordinates.
            [x_sorted, order] = sort(data(i).x);
            
            % Determine what points are on top or bottom surface of airfoil.
            camber_line = NACA_Camber_Line(data(i).x, m, p, t, c, a);
            top_indices = x_sorted' > camber_line;
            bot_indices = ~top_indices;
            
            % Split the order to top indices go first, followed by negative indies in
            % reverse.
            order = [order(top_indices); flip(order(bot_indices))];
            
            % Re-order the elements of the data matrix corresponding to this run.
            data_out(i).IDstr = data(i).IDstr;
            data_out(i).x     = data(i).x(order);
            data_out(i).xnorm = data(i).xnorm(order);
            data_out(i).y     = data(i).y(order);
            data_out(i).p     = data(i).p(order);
            data_out(i).cp    = data(i).cp(order);
            data_out(i).ref   = data(i).ref;
            data_out(i).filename = data(i).filename;
            
        end
        
%     catch
%         
%         delete(hwait);
%         
%     end

    delete(hwait);

end

    
    
    
    
    