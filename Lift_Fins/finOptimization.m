function [optimizedValue,cop,y,cma] = finOptimization(program,constantValue,angle,upperBoundary,nodes)

if program == 1

    % Program '1'

    % It will optimize c_r value with given s and given sweep
    % NOTE [Old Version]: x_t = s; It should be x_t = s * tan (sweep); Since sweep = pi / 4, tan(sweep)=1, so x_t = s
    
    % ------------ INPUT ------------
    
    s = constantValue;
    sweep = angle * (pi/180); % Leave pi/180 as written, just change the sweep angle, it will convert it in radiant.
    
    % Input lower bond and upper bond here, they will be Root Chord
    % Boundaries
    lb = 0; % Lower
    ub = upperBoundary; % Upper. ub must be > s due to geometrial limits. If ub = s, the shape of the fin is a triangle.
    n = nodes; % Number of nodes
    
    sol = zeros(1,n);
    c_t_values = zeros(1,n);
    c_r_values = zeros(1,n);
    x_t = 0;
    x_ott = 0;
    c_r_ott = 0;
    c_t_ott = 0;

    % -------------------------------
    
    if s > 0 && sweep > 0 && n > 0 && lb >= 0 && (ub > 0 && ub > s)
        
        %disp('[Program 1] Optimizing x with c_r...');
        % Defining variables used in this program

        x_t = s * tan(sweep);

        % Function
        function_objective = @(c_r, x_t, c_t) ( x_t / 3 * (c_r + 2 .* c_t)/(c_r + c_t) + 1/6 * (c_r + c_t - (c_r .* c_t)/(c_r + c_t)) );

        c_r_values = linspace(lb, ub, n);

        for i = 1:n
            % I have n c_r values, so I need n c_t values
            c_t_values(i) = c_r_values(i) - x_t;
        end
               
        sol = function_objective(c_r_values, x_t, c_t_values); % It contains all centres of pressure
        
        for i = 1:n
            if sol(i) > x_ott
                x_ott = sol(i);
                c_r_ott = c_r_values(i);
                c_t_ott = c_t_values(i);
            end
        end
        
        optimizedValue = c_r_ott;
        cop = x_ott;
        y = s / 3 * (c_r_ott + 2 * c_t_ott)/(c_r_ott + c_t_ott);
        cma = 2/3 * ( c_r_ott + c_t_ott - (c_r_ott * c_t_ott)/(c_r_ott + c_t_ott) );

        % % ------------ OUTPUT ------------
        % 
        % string_message_c_r_ott = ['Optimized value of c_r is: ', num2str(c_r_ott)];
        % disp(string_message_c_r_ott);
        % 
        % % Calculating centre of pressure with the optimized c_r and given 's'
        % string_message_x_ott = ['Optimized value of x is: ', num2str(x_ott)];
        % disp(string_message_x_ott);
        % 
        % % Calculating spanwise location of the mean aerodynamic chord (y)
        % y_ott = s / 3 * (c_r_ott + 2 * c_t_ott)/(c_r_ott + c_t_ott);
        % string_message_y_ott = ['Optimized value of y is: ', num2str(y_ott)];
        % disp(string_message_y_ott);
        % 
        % % Calculating mean aerodynamic chord length (c_ma)
        % c_ma_ott = 2/3 * ( c_r_ott + c_t_ott - (c_r_ott * c_t_ott)/(c_r_ott + c_t_ott) );
        % string_message_c_ma_ott = ['Optimized value of c_ma is: ', num2str(c_ma_ott)];
        % disp(string_message_c_ma_ott);
        % 
        % % --------------------------------
    else
        if(ub <= s)
            disp('Incorrect initialization of the upper bound [ub]. Code cannot perform.');
        else
            disp('Variables [s, n, sweep, lb, ub] must be positive. Code cannot perform.');
        end
    end

elseif program == 2

    % Program '2'
    
    % This part of the program will optimize the center of pressure
    % location with given 'c_r' and varying 's', which is the Fin Semispan.
    % I will use the same function_objective(c_r, s).

    % ------------ INPUT ------------

    c_r = constantValue;
    sweep = angle * (pi/180); % Leave pi/180 as written, just change the sweep angle, it will convert it in radiant. I know this is redundant, but in order to have a cleaner code I wrote it twice

    lb = 0;
    ub = upperBoundary; % Due to geometrical limits, ub must be < c_r. If it's the same the shape of the fin is a triangle.
    n = nodes;
    
    sol = zeros(1,n);
    s_used = zeros(1,n);
    c_t_used = zeros(1,n);
    x_t_values = zeros(1,n);
    c_t_values = zeros(1,n);
    s_ott = 0;
    x_ott = 0;
    c_t_ott = 0;
    nSolutions = 1;

    % -------------------------------

    if c_r > 0 && sweep > 0 && n > 0 && lb >= 0 && (ub > 0 && ub < c_r)
        % disp('[Program 2] Optimizing x with s...');
        
        % Function
        function_objective = @(c_r, x_t, c_t) ( x_t / 3 * (c_r + 2 .* c_t)/(c_r + c_t) + 1/6 * (c_r + c_t - (c_r .* c_t)/(c_r + c_t)) );
        
        s_values = linspace(lb, ub, n);
        
        for i = 1:n
            % I have n 's' values, so I need n 'x_t' values and n 'c_t'
            % values
            x_t_values(i) = s_values(i) * tan(sweep);
            c_t_values(i) = c_r - x_t_values(i);
        end
        
        for i = 1:n
            x = function_objective(c_r, x_t_values(i), c_t_values(i));
            if x > 0
                sol(nSolutions) = x; % Contains all positive centers of pressure
                s_used(nSolutions) = s_values(i); % It gives the i 's' solution used to have that x_t_values(i) who gave you a positivev 'x'
                c_t_used(nSolutions) = c_t_values(i);
                if nSolutions < n
                    nSolutions = nSolutions + 1;
                end
            end
        end
    
        for i = 1:nSolutions
            if sol(i) > x_ott
                x_ott = sol(i);
                s_ott = s_used(i);
                c_t_ott = c_t_used(i); % It is the one used to calculate x_ott
            end
        end
        
        optimizedValue = s_ott;
        cop = x_ott;
        y = s_ott / 3 * (c_r + 2 * c_t_ott)/(c_r + c_t_ott);
        cma = 2/3 * ( c_r + c_t_ott - (c_r * c_t_ott)/(c_r + c_t_ott) );

        % % ------------ OUTPUT ------------
        % 
        % string_message_c_r = ['Optimized value of s is: ', num2str(s_ott)];
        % disp(string_message_c_r);
        % 
        % % Calculating centre of pressure with the optimized c_r and given 's'
        % string_message_x_ott = ['Optimized value of x is: ', num2str(x_ott)];
        % disp(string_message_x_ott);
        % 
        % % Calculating spanwise location of the mean aerodynamic chord (y)
        % y_ott = s_ott / 3 * (c_r + 2 * c_t_ott)/(c_r + c_t_ott);
        % string_message_y_ott = ['Optimized value of y is: ', num2str(y_ott)];
        % disp(string_message_y_ott);
        % 
        % % Calculating mean aerodynamic chord length (c_ma)
        % c_ma_ott = 2/3 * ( c_r + c_t_ott - (c_r * c_t_ott)/(c_r + c_t_ott) );
        % string_message_c_ma_ott = ['Optimized value of c_ma is: ', num2str(c_ma_ott)];
        % disp(string_message_c_ma_ott);
        % 
        % % --------------------------------   

    else
        if(ub >= c_r)
            disp('Incorrect initialization of the upper bound [ub], it must be minus the [constantValue]. Code cannot perform.');
        else
            disp('Variables [c_r, sweep, lb, ub, n] must be positive. Code cannot perform.')
        end
    end
else
    disp('You chose an incorrect value for variable [program]. It must be 1 or 2.');
end

end