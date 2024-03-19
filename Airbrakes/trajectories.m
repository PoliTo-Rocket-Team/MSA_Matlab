% clear all
% close all
% clc

%% Cd function
T_Cd=readtable("Rocketteam_CD_MD.xlsx");
T_Cd_mat=table2array(T_Cd);
for i=1:11
    Cd=T_Cd_mat(1+20*(i-1):20+20*(i-1),2);
    Mach=T_Cd_mat(1+20*(i-1):20+20*(i-1),3);
    x_Cd=linspace(0,2,100);
    y_Cd{:,i}=@(x) spline(Mach,Cd,x);
    figure(1)
    plot(Mach,Cd,'o')
    fplot(y_Cd,[0 2])
    hold on
end

%% ODE
figure(2)
tspan=[27 8];
y0=[0;3000];
% ft = linspace(0,100,1000);
% f = spline(time,F_t,x);%thrust force from motor

options = odeset('RelTol',1e-3,'Stats','on');

for i=1:11
[t,y{1,i}] = ode45(@(t,y) odefun1(t,y,y_Cd{1,i}), tspan,y0,options);
plot(t,y{1,i}(:,2),'-');
hold on
end
y_mat=cell2mat(y);

function dydt = odefun1(t,y,y_Cd)
M=25.3;
g=9.81;
A=0.0132665;

H=7249;
rho_0=1.225;
rho=@(x) rho_0*exp(-x/H); %x=y(2)

gamma=1.4;
R=287;
T_r=288.04;
K=-0.00649;
T=@(x) T_r+K*x;
c=@(x) sqrt(gamma*R*T(x)); %x=y(2) h


for i=1:11
    Mach=y(1)/c(y(2));
    dydt = zeros(2,1);
    dydt(1)= -g-y_Cd(Mach)*(A/M/2)*rho(y(2))*y(1)^2;
    dydt(2) = y(1); %dh/dt=v
end
end


%% Motor function
% T_m=readtable("Rocketteam_motor.xlsx");
% T_m_mat=table2array(T_m);
% time=[T_m_mat(:,1);linspace(4.6,100,87)'];
% F_t=[T_m_mat(:,2);zeros(87,1)];
% 
% x=linspace(0,100,1000);
% y=spline(time,F_t,x);
% figure(2)
% plot(x,y,'-',time,F_t,'o')