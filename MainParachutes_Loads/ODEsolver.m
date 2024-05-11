close all
clc
clear all

Do = 66/39.37;
Vs = 24.7;
t_f = 8*Do/Vs^0.9;

tspan = [0 t_f];
%y(1)--> v, y(2)--> m_a, y(3)--> F, y(4)--> CdS
y0 = [Vs 0 0 0.1];
[t,y] = ode45(@(t,y) odefun(t,y), tspan, y0);

rho = 1.19934; %Calculation of rho at 270 meters predicted opening of the main chute + 20 meters of ask 't harde
m_p = 0.135;
m_V = 15.801;
p = 0.2;
CdS = 4.8557;
CdS_d = 0.4669;
C_Dr = 0.5;
S_r = 0.01327; %2pi r *l
CdS_r = C_Dr * S_r;
k_g = 1.068*(1-1.465*p-0.25975*p^2+1.2626*p^3);
k_a = ((4/3)*pi*0.6^3*k_g)/((CdS)^(3/2));
W_V = m_V*9.81;
W_p = m_p*9.81;

y_0 = [y(end,1) y(end,2) y(end,3) y(end,4)];
Tspan = [t_f 10];
[t_1,y_1] = ode45(@(t,y) odefun3(t,y), Tspan, y_0);
t = [t;t_1];
y = [y;y_1];

acc = (y(:,3) + 0.5.*rho.*(y(:,1).^2).*(CdS_r+CdS_d) - W_V)./m_V;

figure
plot(t,y(:,1),"Color","red")
title("velocity")
xlabel("t(s)")
ylabel("v(m/s)")
grid on
figure
plot(t,y(:,2),"Color","green")
title("added mass")
xlabel("t(s)")
ylabel("m_a(Kg)")
grid on
figure
plot(t,y(:,3),"color","blue")
title("Force")
xlabel("t(s)")
ylabel("F(N)")
grid on
figure
plot(t,y(:,4),"color","k")
title("CdS")
xlabel("t(s)")
ylabel("CdS(m^2)")
grid on
figure
plot(t,acc,"color","k")
title("Acceleration")
xlabel("t(s)")
ylabel("v'(m/s^2)")
grid on
