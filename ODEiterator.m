close all
tspan = [0 2];
%y(1)--> v, y(2)--> m_a, y(3)--> F, y(4)--> CdS
y0 = [27.379 0 0 0];
t_f_in = linspace(0,2.5,200);
C_Dr_in = linspace(0.2,0.8,7);
p_in = linspace(0,1,11);
figure
% for C_Dr_i = C_Dr_in
%     F_maxes =[];
%     counter=1;
%     for t_f_k = t_f_in
%        [t_ik,y_ik] = ode45(@(t_ik,y_ik) odefun2(t_ik,y_ik, C_Dr_i, 0.2, t_f_k), tspan, y0);
%        F_max=max(y_ik(:,3));
%        F_maxes(counter)=F_max;
%        counter = counter+1;
%     end
%     hold on
%     grid on
%     plot(t_f_in,F_maxes,"-o")
%     test = string(C_Dr_i);
%     text(t_f_in(end),F_maxes(end),test)
% end

% figure
% for p_i = p_in
%     F_maxes =[];
%     counter=1;
%     for t_f_k = t_f_in
%        [t_ik,y_ik] = ode45(@(t_ik,y_ik) odefun2(t_ik,y_ik, 0.5, p_i, t_f_k), tspan, y0);
%        F_max=max(y_ik(:,3));
%        F_maxes(counter)=F_max;
%        counter = counter+1;
%     end
%     hold on
%     grid on
%     plot(t_f_in,F_maxes,"-o")
%     test = string(p_i);
%     text(t_f_in(end),F_maxes(end),test)
% end
F_maxes=[];
counter = 1;
figure
t_fcr=0;
for t_f_k = t_f_in
    [t_ik,y_ik] = ode45(@(t_ik,y_ik) odefun2(t_ik,y_ik, 0.5, 0.2, t_f_k), tspan, y0);
    F_max=max(y_ik(:,3));
    F_maxes(counter)=F_max;
    counter = counter+1;
    if F_max >= 885
        t_f_k=t_fcr
    end
end
    hold on
    grid on
    plot(t_f_in,F_maxes,"-o")
t_fcr
%[t,y] = ode45(@(t,y) odefun(t,y), tspan, y0);


