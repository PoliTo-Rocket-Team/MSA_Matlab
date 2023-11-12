%% LIFT COEFFICNET CALCULATION

clc
clear all

% Characteristic Mach (T=288.25 k, R = 287 J/kgK)
M = 0.0882;
x = linspace(0,pi/10);

% Rocket's radius
r = 52;
% Single fin surface
Sf = 13500; 
% Body surface
Sb = 615124; 

% First method, linearized theory for a thin plate
f1 = @(x) ((2*pi*x)./(sqrt(1-M^2)))*Sf/(2*pi*r^2);

% Second method, lift partition
f2 = @(x) x.*2;

Clf = f1(x);
Cl = f2(Clf);

plot(x*180/pi,Clf,'b');
hold on
plot(x*180/pi,Cl,'r');
hold off

legend({'Cl single fin','Cl rocket'},'Location','northwest');
xlabel('Degree') 
ylabel('Cl') 


%% LIFT COEFFICIENT COMPARISON

clc
clear all

% Characteristic Mach (T=288.25 k, R = 287 J/kgK)
Mc = 0.0882;
% Rocket's radius
r = 52;
% Single fin surface
Sf = 13500; 
% Body surface
Sb = 615124;
%Function values Cl and alpha
[Cl_t,alpha_t] = Simulation;

Sol = zeros(length(Cl_t),4);
i = 1;

for M = 0.3:0.1:0.6
% First method, linearized theory for a thin plate
f1 = @(x) ((2*pi*(x*2*pi)/360)./(sqrt(1-M^2)))*Sf/(2*pi*r^2);

% Second method, lift partition
f2 = @(x) x*2;

Clf = f1(alpha_t);
Cl = f2(Clf);
Sol(:,i) = Cl;
i = i+1;
end

plot(alpha_t,Cl_t,'ob','DisplayName','Cl rocket simulation');
hold on

plot(alpha_t,Sol(:,1),'g','LineWidth',1.5,'DisplayName','Cl rocket approx. M = 0.3');
plot(alpha_t,Sol(:,2),'r','LineWidth',1.5,'DisplayName','Cl rocket approx. M = 0.4');
plot(alpha_t,Sol(:,3),'c','LineWidth',1.5,'DisplayName','Cl rocket approx. M = 0.5');
plot(alpha_t,Sol(:,4),'m','LineWidth',1.5,'DisplayName','Cl rocket approx. M = 0.6');

hold off
legend('Location','northwest')

xlabel('Degree') 
ylabel('Cl') 


%% LIFT COEFFICIENT COMPARISON 2

clc
clear all

% Characteristic Mach (T=288.25 k, R = 287 J/kgK)
Mc = 0.0882;
% Rocket's radius
r = 52;
% Single fin surface
Sf = 13500; 
% Body surface
Sb = 615124;
%Function values Cl and alpha
[Cl_t2,alpha_t2] = Simulation2;

Sol = zeros(length(Cl_t2),4);
i = 1;

for M = 0.50:0.05:0.65
% First method, linearized theory for a thin plate
f1 = @(x) ((2*pi*(x*2*pi)/360)/(sqrt(1-M^2)))*Sf/(2*pi*r^2);

% Second method, lift partition
f2 = @(x) x*2;

Clf = f1(alpha_t2);
Cl = f2(Clf);
Sol(:,i) = Cl;
i = i+1;
end

plot(alpha_t2,Cl_t2,'ob','DisplayName','Cl rocket simulation');
hold on

plot(alpha_t2,Sol(:,1),'g','LineWidth',1.5,'DisplayName','Cl rocket approx. M = 0.5');
plot(alpha_t2,Sol(:,2),'r','LineWidth',1.5,'DisplayName','Cl rocket approx. M = 0.55');
plot(alpha_t2,Sol(:,3),'c','LineWidth',1.5,'DisplayName','Cl rocket approx. M = 0.6');
plot(alpha_t2,Sol(:,4),'m','LineWidth',1.5,'DisplayName','Cl rocket approx. M = 0.65');

hold off
legend('Location','northwest')

xlabel('Degree') 
ylabel('Cl') 


%% BARROWMAN METHOD
clc
clear all

% Dati
C_N_nose = 2;
r = 52;
LF = 119;
CR = 200;
CT = 70;
s = 100;
n = 4;
M = 0.3;
Sf = 13500; 
Sb = 615124;

[Cl_t3,alpha_t3] = Simulation3;

% Barrowman Method
C_N_fins = (4*n*(s/(2*r))^2)/(1+sqrt(1+((2*LF)/(CR+CT))^2));
C_l = @(x)(C_N_nose+C_N_fins)*cos(x*2*pi/360).*((x*2*pi)/360);
C_L = C_l(alpha_t3);

%Barrowman Method with Interference
interference_factor = 1+(r/(s+r));
C_N_fins_interf = C_N_fins * interference_factor;
C_l_interf = @(x)(C_N_nose+C_N_fins_interf)*cos(x*2*pi/360).*((x*2*pi)/360);
C_L_interf = C_l_interf(alpha_t3);

% First Method
f1 = @(x) ((2*pi*(x*2*pi)/360)./(sqrt(1-M^2)))*Sf/(2*pi*r^2);
f2 = @(x) x.*2;
Clf = f1(alpha_t3);
Cl = f2(Clf);

relative_err2 = abs(Cl-Cl_t3)./Cl_t3
relative_err1 = abs(C_L_interf-Cl_t3)./Cl_t3

figure(1)
plot(alpha_t3,relative_err1,'b','DisplayName','Relative Error Barrowman Method')
hold on
plot(alpha_t3,relative_err2,'m','DisplayName','Relative Error First Method')
hold off

legend('Location','northwest')

xlabel('Degree') 
ylabel('Relative Error') 

figure(2)
plot(alpha_t3,Cl_t3,'ob','DisplayName','OpenRocket Simulation');
hold on
plot(alpha_t3,C_L,'r','LineWidth',1.5,'DisplayName','Barrowman Method');
hold on
plot(alpha_t3,C_L_interf,':b','LineWidth',1.5,'DisplayName','Barrowman Method with interference');
hold on
plot(alpha_t3,Cl,'g','LineWidth',1.5,'DisplayName','First Method');
hold off

legend('Location','northwest')

xlabel('Degree') 
ylabel('Cl') 




