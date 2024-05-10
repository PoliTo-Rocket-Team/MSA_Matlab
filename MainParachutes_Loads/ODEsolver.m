close all
tspan = [0 0.7];
%y(1)--> v, y(2)--> m_a, y(3)--> F, y(4)--> CdS
y0 = [27.379 0 0 0];
[t,y] = ode45(@(t,y) odefun(t,y), tspan, y0);


y_0 = [y(end,1) y(end,2) y(end,3) y(end,4)];
Tspan = [0.7 10];
[t_1,y_1] = ode45(@(t,y) odefun3(t,y), Tspan, y_0);
t = [t;t_1];
y = [y;y_1];

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
