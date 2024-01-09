function [optimizedValue,cop,y,cma] = finOptimization(program,initialValue,constantValue,angle)

if program == 1

    % Program '1'

    % It will optimize c_r value with given s and given sweep
    
    % ------------ INPUT ------------
    
    s = constantValue;
    sweep = angle * (pi/180); % [rad]
    
    if s > 0 && sweep > 0
        
        x_t = s * tan(sweep);

        dxdc = @(c_r) ( -(2*c_r^2 - 2*x_t*c_r - x_t^2)/(6*(2*c_r - x_t)^2) );

        [x1,c_r_ott] = fzero(dxdc,initialValue);

        optimizedValue = c_r_ott;
        cop = x1;
        y = s / 3 * (3*c_r_ott - 2*x_t)/(2*c_r_ott-x_t);
        cma = 2/3 * ( 2*c_r_ott - x_t - ( c_r_ott^2 - c_r_ott*x_t )/( 2*c_r_ott - x_t ) );

    else
        
        disp('Variables [s, sweep] must be positive. Code cannot perform.');
        
    end

elseif program == 2

    % Program '2'
    
    % This part of the program will optimize the center of pressure
    % location with given 'c_r' and varying 's', which is the Fin Semispan.

    % ------------ INPUT ------------

    c_r = constantValue;
    sweep = angle * (pi/180); % Leave pi/180 as written, just change the sweep angle, it will convert it in radiant. I know this is redundant, but in order to have a cleaner code I wrote it twice

    % -------------------------------

    if c_r > 0 && sweep > 0
        % disp('[Program 2] Optimizing x with s...');
        
        % Function
        % dxds = @(s) ( tan(sweep) * ( 3*c_r^2-4*c_r*s*tan(sweep)+s^2*(tan(sweep))^2 )/(2*(2*c_r-s*tan(sweep))^2) );
        dxds = @(s) -( tan(sweep) * ( 3*c_r^2-4*c_r*s*tan(sweep)+s^2*(tan(sweep))^2 )/(2*(2*c_r-s*tan(sweep))^2) );
        [x2,s_ott] = fzero(dxds,initialValue);
        
        optimizedValue = s_ott;
        cop = x2;
        y = s_ott / 3 * (3*c_r-2*s_ott*tan(sweep))/(2*c_r-s_ott*tan(sweep));
        cma = 2/3 * ( 2*c_r - s_ott*tan(sweep) - (c_r^2-c_r*s_ott*tan(sweep))/(2*c_r - s_ott * tan(sweep)) );
       
    else
        disp('Variables [c_r, sweep] must be positive. Code cannot perform.')
    end
else
    disp('You chose an incorrect value for variable [program]. It must be 1 or 2.');
end

end