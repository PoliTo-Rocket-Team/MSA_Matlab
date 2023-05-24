=close all
tspan = [0 2];
%y(1)--> v, y(2)--> m_a, y(3)--> F, y(4)--> CdS
y0 = [27.379 0 0 0];
[t,y] = ode45(@(t,y) odefun(t,y), tspan, y0);
figure
plot(t,y(:,1),"Color","red")
title("velocity")
xlabel("t(s)")
ylabel("v(m/s)")
figure
plot(t,y(:,2),"Color","green")
title("added mass")
xlabel("t(s)")
ylabel("m_a(Kg)")
figure
plot(t,y(:,3),"color","blue")
title("Force")
xlabel("t(s)")
ylabel("F(N)")
figure
plot(t,y(:,4),"color","k")
title("CdS")
xlabel("t(s)")
ylabel("CdS(m^2)")
max(y(:,3))
