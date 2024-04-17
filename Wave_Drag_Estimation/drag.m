clear
clc
close all

%% DEFINITION AND SETUP OF VARIABLES
mach = 2.0;
air_density = 1.225;
velocity = mach * sqrt(1.4 * 287 * 288);
[root_chord,tip_chord,span,thickness,sweep_angle,x_le,x_te_t,...
    x_te_r,c_1,c_2] = setup_variables; % look up the fun to change values

%% MULTIPOLE STRENGTH FUNCTIONS FOR THE FINS
beta = sqrt(mach.^2 - 1);

% theta = linspace(0, 2 * pi,1000);
% elemental_area_0 = zeros(1,length(theta));
% elemental_area_n = zeros(1,length(theta));
% for i = 1:length(theta)
%     y = linspace(0,span,10);
%     func_val_0 = zeros(1,length(y));
%     func_val_n = zeros(1,length(y));
%     for j = 1:length(y)
%         func_val_0(j) = thickness_wing(x + beta * y(j) * ...
%             cos(theta(i)),y(j));
%         func_val_n(j) = thickness_wing(x + beta * y(j) * ...
%             cos(theta(i)),y(j)) * cos(n * theta(i));
%     end
%     elemental_area_n(i) = trapz(y,func_val_n);
%     elemental_area_0(i) = trapz(y,func_val_0);
% end
% F_0_w_n = (1 / (2 * pi)) * trapz(theta,elemental_area_0);
% F_n_w_n = (1 / pi) * trapz(theta,elemental_area_n)

% x_2_bound = @(x_1) x_1;
% drag_0 = -((0.5 * air_density * velocity^2) / (2 * pi)) * ...
%     (integral2(@(x_1,x_2) final_integrand_0(x_1,x_2,beta), 0, 0.355, ...
%     x_2_bound, 0.355) + integral2(@(x_1,x_2) ...
%     final_integrand_0(x_1,x_2,beta), 0, 0.355, 0, x_2_bound))
% x_1 = linspace(0,0.355,300);
% x_2 = linspace(0,0.355,100);
% f1 = d2F_0_w(x_1,beta);
% f2 = d2F_n_w(x_1,beta,2);
% plot(x_1,f1)
% figure
% plot(x_1,f2)
% x_interp = linspace(0,root_chord,300);
% y_interp = d2F_0_w(x_interp,beta);
% f_0_w = @(x) interp1(x_interp,y_interp,x);
% upper_bound = @(x) x - 3e-1;
% lower_bound = @(x) x + 1e-1;
% x_integ = linspace(0,root_chord,5000);
% y_integ = linspace(0,root_chord,5000);
% [X,Y] = meshgrid(x_integ,y_integ);
% F_0 = final_integrand_0(X,Y,beta,300);
% drag_out = trapz(y_integ,trapz(x_integ,F_0,2));
% i = 2;
% for n = 1:i
%     F_N = final_integrand_n(X,Y,beta,2*n,300);
%     drag_out = drag_out + 0.5*trapz(y_integ,trapz(x_integ,F_N,2));
% end
% drag_out = drag_out * (- 0.5 * air_density * velocity^2) / (2 * pi)
theta = linspace(0,2*pi,300);
y = final_integrand(theta,beta) * ...
    -(0.5 * air_density * velocity^2) / (2 * pi);
drag_w = trapz(theta,y) / (2 * pi);
c_d_w = drag_w / (0.5 * air_density * velocity^2 * 0.013273)

%% setup of the variables needed for various functions
function [root_chord,tip_chord,span,thickness,sweep_angle,x_le,x_te_t,...
    x_te_r,c_1,c_2] = setup_variables
    root_chord = 0.355; % in m
    tip_chord = 0.194; % in m
    span = 0.085; % in m
    thickness = 0.008; % in m
    sweep_angle = 60; % in deg
    root_leading_edge_ratio = 0.2;
    root_trailing_edge_ratio = 0.2;
    sweep_angle = deg2rad(sweep_angle);
    x_le = span * tan(sweep_angle);
    c_1 = root_chord * root_leading_edge_ratio;
    c_2 = root_chord * root_trailing_edge_ratio;
    x_te_t = x_le + tip_chord - c_2;
    x_te_r = root_chord - c_2;
end

%% thickness function definition
function t = thickness_wing(x,y)
    [root_chord,~,span,thickness,~,x_le,x_te_t,x_te_r,...
        c_1,c_2] = setup_variables;
    t = zeros(1,length(y)); %thickness

% recap of functions of the fin, the function itself is made to take
% arrays in both the x and y coordinates, with values with the same index
% from both arrays corresponding to a unique point in the (x,y) plane

% leading_edge_thickness = @(x,y) (x - y * x_le / span) * thickness / c_1;
% min_x_leading = @(y) y * x_le / span;
% max_x_leading = @(y) y * x_le / span + c_1;
% planform_thickness = @(x,y) thickness;
% min_x_planform = @(y) y * x_le / span + c_1;
% max_x_planform = @(y) y * (x_te_t - x_te_r) / span + x_te_r;
% trailing_edge_thickness = @(x,y) (- x + y * (x_te_t - x_te_r) / span ...
%     + root_chord ) * thickness / c_2;
% min_x_trailing = @(y) y * (x_te_t - x_te_r) / span + x_te_r;
% max_x_trailing = @(y) y * (x_te_t - x_te_r) / span + root_chord;

    for i = 1:length(y)
        if y(i) >= 0 && y(i) <= span
            if y(i) * x_le / span <= x(i) && y(i) * x_le / span + c_1 ...
                    >= x(i)
                t(i) = (x(i) - (y(i) * x_le / span)) * thickness / c_1;
            elseif y(i) * x_le / span + c_1 <= x(i) && y(i) * (x_te_t - ...
                    x_te_r) / span + x_te_r >= x(i)
                t(i) = thickness;
            elseif y(i) * (x_te_t - x_te_r) / span + x_te_r <= x(i) && ...
                    y(i) * (x_te_t - x_te_r) / span + root_chord >= x(i)
                t(i) = (- x(i) + (y(i) * (x_te_t - x_te_r) / span + ...
                    root_chord) ) * thickness / c_2;
            else
                t(i) = 0;
            end
        end
    end
end

%% definition of integrand function for multipole strength calculation
function sw = elemental_area_integrand_0(x,theta,beta)
    [~,~,span] = setup_variables;
    elemental_area = @(x,theta_p,beta) integral(@(y) thickness_wing(x + ...
        beta * y * cos(theta_p), y), 0, span);
    sw = zeros(1,length(theta));
    for i = 1:length(theta)
        sw(i) = elemental_area(x,theta(i),beta);
    end
end


function d2f = d2_elemental_area(x,theta,beta)
    [~,~,span] = setup_variables;
    elemental_area = @(x,theta_p,beta) integral(@(y) thickness_wing(x + ...
        beta * y * cos(theta_p), y), 0, span);
    f = zeros(1,length(x));
    for i = 1:length(x)
        f(i) = elemental_area(x(i),theta,beta);
    end
    df = gradient(f,x);
    d2f = gradient(df,x);
end


function sw = elemental_area_integrand_n(x,theta,beta,n)
    [~,~,span] = setup_variables;
    elemental_area = @(x,theta_p,beta) integral(@(y) thickness_wing(x + ...
        beta * y * cos(theta_p), y), 0, span) * cos(n * theta_p);
    sw = zeros(1,length(theta));
    for i = 1:length(theta)
        sw(i) = elemental_area(x,theta(i),beta);
    end
end


% function d2f = d2F_0_w_g(x, beta)
%     F_0_w_f = @(x,beta) (1 / (2 * pi)) * integral(@(theta) ...
%     elemental_area_integrand_0(x,theta,beta), 0, 2 * pi);
%     f = zeros(1,length(x));
%     for i = 1:length(x)
%         f(i) = F_0_w_f(x(i),beta);
%     end
%     d1f = gradient(f)./gradient(x);
%     d2f = gradient(d1f)./gradient(x);
% end


function d2f = d2F_0_w(x,beta)
    F_0_w_f = @(x,beta) (1 / (2 * pi)) * integral(@(theta) ...
    elemental_area_integrand_0(x,theta,beta), 0, 2 * pi);
    f = zeros(1,length(x));
    for i = 1:length(x)
        f(i) = F_0_w_f(x(i),beta);
    end
    d2f = zeros(1,length(x));
    for i = 2:length(x)-1
        d2f(i) = ((f(i + 1) - f(i)) / (x(i + 1)-x(i)) - ...
            (f(i) - f(i - 1)) / (x(i) - x(i - 1))) / ...
            (0.5 * (x(i + 1) - x(i - 1)));
    end
    d2f(1) = d2f(2);
    d2f(end) = d2f(end-1);
end


function d2f = d2F_n_w(x,beta,n)
    F_0_w_f = @(x,beta) (1 / pi) * integral(@(theta) ...
    elemental_area_integrand_n(x,theta,beta,n), 0, 2 * pi);
    f = zeros(1,length(x));
    for i = 1:length(x)
        f(i) = F_0_w_f(x(i),beta);
    end
    d2f = zeros(1,length(x));
    for i = 2:length(x)-1
        d2f(i) = ((f(i + 1) - f(i)) / (x(i + 1)-x(i)) - ...
            (f(i) - f(i - 1)) / (x(i) - x(i - 1))) / ...
            (0.5 * (x(i + 1) - x(i - 1)));
    end
    d2f(1) = d2f(2);
    d2f(end) = d2f(end-1);
% d2f(i) = ((f(i+1)-f(i))/(x(i+1)-x(i)) - (f(i)-f(i-1))/(x(i)-x(i-1)))/...
% (0.5*(x(i+1)-x(i-1)));
end


function res = final_integrand_0_punctual(x_1,x_2,beta)
    f1 = d2F_0_w(x_1(1,:), beta);
    f2 = d2F_0_w(x_2(:,1), beta);
    [F1,F2] = meshgrid(f1,f2);
    res = zeros(size(x_1));
    for i = 1:length(F1(1,:))
        for j = 1:length(F2(:,1))
            if x_1(j,i) == x_2(j,i)
                res(j,i) = 0;
            else
                res(j,i) = F1(j,i) * F2(j,i) * ...
                    log(abs(x_1(j,i) - x_2(j,i)));
            end
        end
    end
end


function d2f = d2S_w(x,theta,beta)
    [~,~,span] = setup_variables;
    elemental_area = @(x,theta,beta) integral(@(y) thickness_wing(x + ...
        beta * y * cos(theta), y), 0, span);
    f = zeros(1,length(x));
    for i = 1:length(x)
        f(i) = elemental_area(x(i),theta,beta);
    end
    d2f = zeros(1,length(x));
    for i = 2:length(x)-1
        d2f(i) = ((f(i + 1) - f(i)) / (x(i + 1)-x(i)) - ...
            (f(i) - f(i - 1)) / (x(i) - x(i - 1))) / ...
            (0.5 * (x(i + 1) - x(i - 1)));
    end
    d2f(1) = d2f(2);
    d2f(end) = d2f(end-1);
end


function int = integrand_d(x_1,x_2,theta,beta)
    f1 = d2S_w(x_1(1,:),theta,beta);
    f2 = d2S_w(x_2(:,1),theta,beta);
    [F1,F2] = meshgrid(f1,f2);
    int = zeros(size(x_1));
    for i = 1:length(F1(1,:))
        for j = 1:length(F2(:,1))
            if x_1(j,i) == x_2(j,i)
                int(j,i) = 0;
            else
                int(j,i) = F1(j,i) .* F2(j,i) .* ...
                    log(abs(x_1(j,i) - x_2(j,i)));
            end
        end
    end
end


function res = final_integrand_0(x_1,x_2,beta,n_int)
    [root_chord] = setup_variables;
    x_interp = linspace(0,root_chord,n_int);
    y_interp = d2F_0_w(x_interp,beta);
    f_0_w = @(x) interp1(x_interp,y_interp,x);
    integrand = @(x_1,x_2) f_0_w(x_1).*f_0_w(x_2).*log(abs(x_1 - x_2));
    res = integrand(x_1,x_2);
    for i = 1:length(res(:,1))
        for j = 1:length(res(1,:))
            if isnan(res(i,j)) || res(i,j) < -100
                res(i,j) = -100;
            end
        end
    end
end


function res = final_integrand_n(x_1,x_2,beta,n,n_int)
    [root_chord] = setup_variables;
    x_interp = linspace(0,root_chord,n_int);
    y_interp = d2F_n_w(x_interp,beta,n);
    f_n_w = @(x) interp1(x_interp,y_interp,x);
    integrand = @(x_1,x_2) f_n_w(x_1).*f_n_w(x_2).*log(abs(x_1 - x_2));
    res = integrand(x_1,x_2);
    for i = 1:length(res(:,1))
        for j = 1:length(res(1,:))
            if isnan(res(i,j)) || res(i,j) < -100
                res(i,j) = -100;
            end
        end
    end
end


function res = another_integrand(x_1,x_2,theta,beta,n_int)
    [root_chord] = setup_variables;
    x_interp = linspace(0,root_chord,n_int);
    y_interp = d2_elemental_area(x_interp,theta,beta);
    d2f_int = @(x) interp1(x_interp,y_interp,x);
    integrand = @(x_1,x_2) d2f_int(x_1).*d2f_int(x_2).*log(abs(x_1 - x_2));
    res = integrand(x_1,x_2);
    for i = 1:length(res(:,1))
        for j = 1:length(res(1,:))
            if isnan(res(i,j)) || res(i,j) < -20
                res(i,j) = -20;
            end
        end
    end
end


function res = final_integrand(theta,beta)
    [root_chord] = setup_variables;
    res = zeros(1,length(theta));
    x_integ = linspace(0,root_chord,1000);
    [X_1,X_2] = meshgrid(x_integ,x_integ);
    for i = 1:length(theta)
        F = another_integrand(X_1,X_2,theta(i),beta,300);
        res(i) = trapz(x_integ,trapz(x_integ,F,2));
    end
end